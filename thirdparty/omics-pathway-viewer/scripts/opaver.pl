#!/usr/bin/env perl
# parse_prot_tran_data.pl
# ver 0.1
# 2014/09/13
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.
#
# 2017/07/21
# - add support for metabolic expression data
# 2015/05/23
# - add codes to take average proteomic data if genes has multiple expression level
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use LWP::UserAgent;
use JSON;
use strict;

$|=1;
my %opt;
my $res=GetOptions( \%opt,
    'prot|a=s',
    'prot_pv_cutoff|ap=s',
    'tran|t=s',
    'tran_pv_cutoff|tp=s',
    'meta|m=s',
    'meta_pv_cutoff|mp=s',
    'outdir|o=s',
    'prefix|p=s',
    'koinfojson|k=s',
    'pwyinfojson|y=s',
    'help|?') || &usage();

if ( $opt{help} || ! (-e  $opt{tran} || -e  $opt{prot} || -e  $opt{meta}) ) { &usage(); }

# Create a user agent
my $ua = LWP::UserAgent->new();
$ua->env_proxy;

#init
my $prefix  = defined $opt{prefix} ? $opt{prefix} : "output";
my $outdir  = defined $opt{outdir} ? $opt{outdir} : "./$prefix";
my $prot_p  = defined $opt{prot_pv_cutoff} ? $opt{prot_pv_cutoff} : 0.05;
my $tran_p  = defined $opt{tran_pv_cutoff} ? $opt{tran_pv_cutoff} : 0.05;
my $meta_p  = defined $opt{meta_pv_cutoff} ? $opt{meta_pv_cutoff} : 0.05;
my $kojson  = defined $opt{koinfojson} ? $opt{koinfojson} : "$outdir/info_ko.json";
my $pwyjson = defined $opt{pwyinfojson} ? $opt{pwyinfojson} : "$outdir/info_pathway.json";
my ($koinfo,$pwyinfo,$map,$pathway,$gene_exp,$geneid2ko);

#init directory
`mkdir -p $outdir`;
die "[ERROR] Can't create output directory $outdir: $!\n" unless -d $outdir;

#init KEGG info
$koinfo = &retrieveKoInfoJson( $kojson ) if -s $kojson;

#init pathway info
if( -s $pwyjson ){
	$pwyinfo = &retrievePathwayInfoJson($pwyjson);
}
if( ! defined $pwyinfo->{"00010"} ){
	$pwyinfo =  &builtPathwayInfoJson();
	&writeJSON($pwyinfo, $pwyjson);
}

my $time = time;

print STDERR "[INFO] Parsing TRAN files...\n";
my ($tran, $t_pos, $t_neg);
if( -e $opt{tran} ){
	($tran, $t_pos, $t_neg) = &parseExpressionFile($opt{tran}, $tran_p, "tran"); 
	print STDERR "[INFO] Done.\n";
}
else{
	print STDERR "[INFO] No input data for transcriptomics. Skipped.\n";
}

print STDERR "[INFO] Parsing PROT files...\n";
my ($prot, $p_pos, $p_neg);
if( -e $opt{prot} ){
	($prot, $p_pos, $p_neg) = &parseExpressionFile($opt{prot}, $prot_p, "prot");
	print STDERR "[INFO] Done.\n";
}
else{
	print STDERR "[INFO] No input data for proteomics. Skipped.\n";
}


print STDERR "[INFO] Parsing METAB files...\n";
my ($meta, $m_pos, $m_neg);
if( -e $opt{meta} ){
	($meta, $m_pos, $m_neg) = &parseExpressionFile($opt{meta}, $meta_p, "meta");
	print STDERR "[INFO] Done.\n";
}
else{
	print STDERR "[INFO] No input data for metabolism. Skipped.\n";
}


foreach my $ko ( keys %$tran ){
	my @mapid = split /, /, $koinfo->{$ko}->{pathway};
	foreach my $mapid (@mapid){
		$map->{$mapid}->{$ko}->{tran} = $tran->{$ko};
		$map->{$mapid}->{$ko}->{info} = $koinfo->{$ko};
	}
}

foreach my $ko ( keys %$prot ){
	my @mapid = split /, /, $koinfo->{$ko}->{pathway};
	foreach my $mapid (@mapid){
		$map->{$mapid}->{$ko}->{prot} = $prot->{$ko};
		$map->{$mapid}->{$ko}->{info} = $koinfo->{$ko};
	}
}

foreach my $ko ( keys %$meta ){
	my @mapid = split /, /, $koinfo->{$ko}->{pathway};
	foreach my $mapid (@mapid){
		$map->{$mapid}->{$ko}->{cpd} = $meta->{$ko};
		$map->{$mapid}->{$ko}->{info} = $koinfo->{$ko};
	}
}

print STDERR "Writing tran/prot expression data...";
open GOUT, ">$outdir/exp_gene.txt" or die "[ERROR] Can't write expression file: $!\n";
open MOUT, ">$outdir/exp_cpd.txt" or die "[ERROR] Can't write expression file: $!\n";
print GOUT "KO\tPROT_GENE\tTRAN_GENE\tPROT\tTRAN\n";
print MOUT "KO\tCOMPOUND\tFC\n";

foreach my $ko ( keys %$gene_exp ){
	next unless $ko;
	if( defined $gene_exp->{$ko}->{meta} ){
		my $m  = $gene_exp->{$ko}->{meta}->{val};
		my $mn = $gene_exp->{$ko}->{meta}->{name};
		if( $m ne "" ){
			printf MOUT "%s\t%s\t%s\n",
				$ko,
				defined $mn ? $mn : "",
				defined $m  ? $m : "";
		}
	}
	else{
		my $p  = $gene_exp->{$ko}->{prot}->{val};
		my $t  = $gene_exp->{$ko}->{tran}->{val};
		my $pg = $gene_exp->{$ko}->{prot}->{name};
		my $tg = $gene_exp->{$ko}->{tran}->{name};
		if( $p ne "" || $t ne "" ){
			printf GOUT "%s\t%s\t%s\t%s\t%s\n",
				$ko,
				defined $pg ? $pg : "",
				defined $tg ? $tg : "",
				#defined $p  ? ($p>0?$p/$p_pos:$p/$p_neg) : "",
				#defined $t  ? ($t>0?$t/$t_pos:$t/$t_neg) : "";
				defined $p  ? $p : "",
				defined $t  ? $t : "";
		}
	}
}
close GOUT;
close MOUT;
print STDERR "Done\n";

print STDERR "Writing pathway list...";
open OUT, ">$outdir/exp_pathway.txt" or die "[ERROR] Can't write pathway file: $!\n";
foreach my $mapid ( sort {$pathway->{$b}<=>$pathway->{$a}} keys %$pathway ){
	print OUT "$mapid\t$pathway->{$mapid}\t$pwyinfo->{$mapid}\n";
}
close OUT;
print STDERR "Done\n";

print STDERR "Writing expression JSON files for KEGG maps...";
foreach my $mapid ( keys %$map ){
    &writeJSON($map->{$mapid}, "$outdir/ko$mapid.exp.json");
}
print STDERR "Done\n";

print STDERR "Writing KO JSON files...";
&writeJSON($koinfo, $kojson);
print STDERR "Done\n";

####################################################################################

sub parseExpressionFile {
	my ($file, $pval_cutoff, $type) = @_;
	my $data;
	my $max_pos=0;
	my $max_neg=0;

	open EXP, $file or die "[ERROR] Can't open expression data: $!\n";
	foreach(<EXP>){
		chomp;
		#skip header
		next if /^NAME/i;
		next if /Pval/;
		next if /^#/;

		my ($gene,$kostr,$logfc,$pvalue) = split /\t/, $_;
		$pvalue = 0 if $pvalue eq "";
		next if !$kostr;
		
		$kostr =~ s/["',;]/ /g;
		my @temp = split /\s+/, $kostr;

		foreach my $ko ( @temp ){
			next if $ko !~ /^(K|C)/;
			$geneid2ko->{$gene} = $ko;
			next if $pvalue > $pval_cutoff;
			
			if($logfc > 0){
				$max_pos = $logfc if $logfc > $max_pos;
			}
			else{
				$max_neg = abs($logfc) if abs($logfc) > $max_neg;
			}

			$data->{$ko}->{$gene}->{logfc}=$logfc;
			$data->{$ko}->{$gene}->{pvalue}=$pvalue;

			&getKeggInfo($ko);
			
			if( abs($logfc) > abs($gene_exp->{$ko}->{$type}->{val}) ){
				$gene_exp->{$ko}->{$type}->{val} = $logfc;
				$gene_exp->{$ko}->{$type}->{name}= $gene;
			}
			print STDERR "$gene\t$ko\t$logfc\n";

			die "[ERROR] Conflict geneId <=> ko mapping: $gene => $geneid2ko->{$gene} and $ko!\n" if defined $geneid2ko->{$gene} && $geneid2ko->{$gene} ne $ko; 
		}
	}
	close EXP;
	return ($data,$max_pos,$max_neg);
}

sub log2 {
	my $v = shift;
	return (log($v)/log(2));
}

sub writeJSON {
    my ($map, $outfile) = @_;
    open OUT, ">$outfile" or die "[ERROR] Can't write JSON file: $!\n";
    my $json_text = to_json($map, {utf8 => 1, pretty => 1});
    print OUT $json_text;
    close OUT;
}

sub retrieveKoInfoJson {
	my $jsonfile = shift;
	open JSON, $jsonfile or die "[ERROR] Can't read KO info JSON file: $!\n";
	local $/ = undef;
	my $json_text = from_json(<JSON>);
	close JSON;
	return $json_text;
}

sub retrievePathwayInfoJson {
	my $jsonfile = shift;
	open JSON, $jsonfile or die "[ERROR] Can't read Pathway info JSON file: $!\n";
	local $/ = undef;
	my $json_text = from_json(<JSON>);
	close JSON;
	return $json_text;
}

sub builtPathwayInfoJson {
	# URL for service
	my $url = "http://rest.kegg.jp/list/pathway";
	my $response = $ua->post("$url");
	my $pwy;
	
	# Check for HTTP error codes
	print STDERR "ERROR: http://rest.kegg.jp/list/pathway\n" unless ($response->is_success); 
	  
	# Output the entry
	my $content = $response->content();
	my @lines = split /\n/, $content;
	
	foreach my $line (@lines){
		my ($id,$name) = $line =~ /path:map(\d+)\t(.*)$/;
		$pwy->{$id} = $name;
	}
	return $pwy;
}

sub getKeggInfo {
	my $ko = shift;
	my @p;

	if( !defined $koinfo->{$ko} ){
		# URL for service
	    my $info;
	    my $url = "http://rest.kegg.jp/get";
		my $response = $ua->post("$url/$ko");

		# Check for HTTP error codes
		unless ($response->is_success){
			print STDERR "No $ko info found.\n";
			return;
		}
		  
		# Output the entry
		my $content = $response->content();
		my @lines = split /\n/, $content;
	
		my $pwyflag=0;
		foreach my $line (@lines){
			last if $line =~ /^\w+/ && $pwyflag;
			$pwyflag = 1 if $line =~ /PATHWAY\s+/;
			if( $pwyflag ){
				# koXXXXXX for KO, mapXXXXXX for compound
				if( $line =~ /\s+ko(\d+)/ || $line =~ /\s+map(\d+)/ ){
					push @p, $1; 
					$pathway->{$1} = 1;
				}
			}
	        else{
	            $info->{name}       = $1 if $line =~ /NAME\s+(.*);?$/;
	            $info->{formula}    = $1 if $line =~ /FORMULA\s+(.*)$/;
	            $info->{definition} = $1 if $line =~ /COMMENT\s+(.*)$/; # for compound record, e.g. http://rest.kegg.jp/get/C01290
	            $info->{definition} = $1 if $line =~ /DEFINITION\s+(.*)$/;
	        }
		}
	    $info->{pathway} = join ", ", @p;
		$koinfo->{$ko} = $info;
	}

	my $info = $koinfo->{$ko};
	my @mapids = split /, /, $info->{pathway};
	foreach my $mapid ( @mapids ){
		$pathway->{$mapid} ||= 0;
		$pathway->{$mapid}++;
	}
}

sub timeInterval{
    my $now = shift;
    $now = time - $now;
    return sprintf "%02d:%02d:%02d", int($now / 3600), int(($now % 3600) / 60), int($now % 60);
}

sub usage {
print <<__END__;
OPaver (Omics Pathway Viewer) is a web-based viewer to provide 
an integrated and interactive visualization of omics data in KEGG maps.

Version: v0.3.1

$0 [OPTIONS] -a <PROT> -t <TRAN> -m <META>
	--prot|a              proteomic expression data in TSV format
	--tran|t              transcriptomic expression data in TSV format
	--meta|m              metabolic expression data in TSV format

[OPTIONS]
	--prot_pv_cutoff|ap   pvalue cutoff, default is 0.05
	--tran_pv_cutoff|tp   pvalue cutoff, default is 0.05
	--meta_pv_cutoff|mp   pvalue cutoff, default is 0.05
	--outdir|o
	--prefix|p
	--koinfojson|k
	--pwyinfojson|y

__END__
exit();
}
