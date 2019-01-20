package cocosv;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::Perl;


sub feature_types{
	return['Transcript']
}

sub get_header_info {
    my $self = shift;
    return {
	EJUNCT_POS => 'Position of Last Exon-Exon Junction',
	REF_SEQ => 'Coding sequence of reference transcript',
	REF_SEQ_PEPTIDE => 'Peptide sequence of reference transcript',
	VAR_SEQ => 'Coding sequence of variant transcript',
	VAR_SEQ_PEPTIDE => 'Peptide sequence of variant transcript'
    };
}


sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);
    return $self;
}


sub translate_seq {
	
	## Add checks for incomplete annotations at 5'
	my $seq = shift;
	return "" if $seq eq "";

	my %codon_tbl = (
		"TTT" => "F",
		"TTC" => "F",
		"TTA" => "L",
		"TTG" => "L",
		"CTT" => "L",
		"CTC" => "L",
		"CTA" => "L",
		"CTG" => "L",
		"ATT" => "I",
		"ATC" => "I",
		"ATA" => "I",
		"ATG" => "M",
		"GTT" => "V",
		"GTC" => "V", 
		"GTA" => "V", 
		"GTG" => "V",
		"TCT" => "S",
		"TCC" => "S",
		"TCA" => "S", 
		"TCG" => "S", 
		"CCT" => "P",
		"CCC" => "P", 
		"CCA" => "P", 
		"CCG" => "P", 
		"ACT" => "T", 
		"ACC" => "T", 
		"ACA" => "T", 
		"ACG" => "T", 
		"GCT" => "A", 
		"GCC" => "A", 
		"GCA" => "A", 
		"GCG" => "A",
		"TAT" => "Y", 
		"TAC" => "Y", 
		"TAA" => "STOP", 
		"TAG" => "STOP", 
		"CAT" => "H", 
		"CAC" => "H", 
		"CAA" => "Q", 
		"CAG" => "Q", 
		"AAT" => "N", 
		"AAC" => "N", 
		"AAA" => "K", 
		"AAG" => "K", 
		"GAT" => "D", 
		"GAC" => "D", 
		"GAA" => "E", 
		"GAG" => "E", 
		"TGT" => "C",
		"TGC" => "C", 
		"TGA" => "STOP", 
		"TGG" => "W", 
		"CGT" => "R", 
		"CGC" => "R", 
		"CGA" => "R", 
		"CGG" => "R", 
		"AGT" => "S", 
		"AGC" => "S", 
		"AGA" => "R", 
		"AGG" => "R", 
		"GGT" => "G", 
		"GGC" => "G", 
		"GGA" => "G", 
		"GGG" => "G"
	);

	if($seq =~ m/ATG/){
		
		#print "\n *** Initial Sequence *** \n";
		#print $seq, "\n";
		my $local_seq = substr($seq, $-[0]);
		my $tr = "";
		my $str_idx = 0;
		while($str_idx < length($local_seq) - 1){
			my $local_trs = $codon_tbl{substr($local_seq, $str_idx, 3)};
			last if $local_trs eq 'STOP';
			$tr.=$local_trs;
			$str_idx+=3;
		}
		return $tr;
	}

	return "";	
}

sub get_seqs {
	my ($tva) = @_;
	my $trans = $tva->transcript;
	my $trans_seq = $trans->spliced_seq;
	my @coding_coord = ($trans->cdna_coding_start(),$trans->cdna_coding_end());

	if(defined $coding_coord[0] && defined $coding_coord[1]){
		## Get All Exon Positions ##
		my @all_exon = @{$trans->get_all_Exons()};
		my @all_exon_end_pos = map{$_->cdna_end($trans)}(@all_exon);
		my @all_exon_start_pos = map{$_->cdna_start($trans)}(@all_exon);
		
		## Get the last exon with coding region
		my $idx = 0;
		while($idx < (@all_exon_end_pos-1) && $all_exon_end_pos[$idx] < $coding_coord[1]){
			$idx+=1;
		}
		my $ejunct_pos = $all_exon_end_pos[$idx-1]; 
		
		my $ref_coding_seq = substr($trans_seq, $coding_coord[0], $coding_coord[1]);
		print "\n\n** Reference Coding Sequence **\n\n";
		print $ref_coding_seq, "\n";
		my $allele_str = $tva->allele_string;
		my $var_pos_start = $tva->transcript_variation->cds_start();
		my $var_pos_end = $tva->transcript_variation->cds_end();
		
		## Modify Sequence
		my $var_coding_seq = "";
		if(defined $var_pos_start && defined $var_pos_end){
			my $var_allele = (split /\//, $allele_str)[1];
			$var_allele =~ s/-//;
			$var_coding_seq = substr($ref_coding_seq, 0, $var_pos_start - 1).$var_allele.substr($ref_coding_seq, $var_pos_end-1);
			print "\n\n** Variant Coding Sequence **\n\n";
			print $var_coding_seq, "\n";
		}
		
		## Translate Sequences ##
		my $ref_peptide_seq = translate_seq($ref_coding_seq);
		my $var_peptide_seq = translate_seq($var_coding_seq);
		return($ref_coding_seq,$var_coding_seq,$ref_peptide_seq,$var_peptide_seq,$ejunct_pos)
	}
	return "";
}


sub run {
	my ($self, $tva, $vep_line) = @_;
	my ($ref_coding_seq,$var_coding_seq,$ref_peptide_seq,$var_peptide_seq,$ejunct_pos) = get_seqs($tva);

	return {
		EJUNCT_POS => $ejunct_pos,
		REF_SEQ => $ref_coding_seq,
		REF_SEQ_PEPTIDE => $ref_peptide_seq,
		VAR_SEQ => $var_coding_seq,
		VAR_SEQ_PEPTIDE => $var_peptide_seq
	};
}

1;

