SCOPe: Structural Classification of Proteins — extended. Release 2.08 (updated 2023-01-06, stable release September 2021)
References: Fox NK, Brenner SE, Chandonia JM. 2014. Nucleic Acids Research 42:D304-309. doi: 10.1093/nar/gkt1240.
Chandonia JM, Guan L, Lin S, Yu C, Fox NK, Brenner SE. 2022. Nucleic Acids Research 50:D553–559. doi: 10.1093/nar/gkab1054. (citing information)

からのデータを
ESM2-650M 
Rives, Alexander, et al. "Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences." Proceedings of the National Academy of Sciences 118.15 (2021): e2016239118.
にかけて特徴量を抽出した。

TMALIGN のアラインメント
>d6iyia_.pdb
SLTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSSA-QVHTHGKKIIDALTTAVNHID-DI-TG-ALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR------
>d7diha_.pdb
MLSEETIRVIKSTVPLLKEHGTEITARMFELLFSKYPKTKELFAGA--------SE--EQPKKLANAIIAYATYIDRLEELDNAISTIARSHV-RRNVKPEHYPLVKECLLQAIEEVLN-P-GEEVLKAWEEAYDFLAKTLITLEKKLYSQP

>d6iyia_.pdb
SLTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSS--A-QVHTHGKKIIDALTTAVNHIDDITGALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR
>d4f6ia_.pdb
---------------MV--NWAAVVDDFYQELFKAHPEYQNKFGFKGVALGSLKGNAAYKTQAGKTVDYINAAIGGS--A--DAAGLASRHKG-RNVGSAEFHNAKACLAKACSAHG---A---PD--LGWAIDDILSH-L----

そもそも < 0.5 
>d6iyia_.pdb
S------------------------------LTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSSAQVHTHGKKIIDALTTAVNHIDDI-TGALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR
>d6lumb2.pdb
-MEPFFDAYRAVKPFLVTSGNPPTKERIQSPTDRARYD-DT--TK-CILC-ACCTTSC--PVYW--S--E--G--S---------Y--FG--PAAIVNAHRFIFDS-RDEAAAERLDILN--E-VDGVWRCR--TT-FNCTEACP--RG---IQ-VTQ-AIQEVKRALMFA--



>seq1
SLTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSSAQVHTHGKKIIDALTTAVNHIDD---ITGALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR------
>d7diha_.res.gz
MLSEETIRVIKSTVPLLKEHGTEITARMFELLFSKYPKTKELFAGAS---------EEQPKKLANAIIAYATYIDRLEELDNAISTIARSHV-RRNVKPEHYPLVKECLLQAIEEVLN--PGEEVLKAWEEAYDFLAKTLITLEKKLYSQP

>seq1
SLTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADL---SVSSAQVHTHGKKIIDALTTAVNHIDDITGALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR
>d4f6ia_.res.gz
M-----------------VNWAAVVDDFYQELFKAHPEYQNKFGFKGVALGSLKGNAAYKTQAGKTVDYINAAIG----GSADAAGLASRHK-GRNVGSAEFHNAKACLAKACSAHGAPDLGWAIDDI-------------LSHL

>seq1
SLTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSSAQVHTHGKKIIDALTTAVNHIDDITGALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR
>d6lumb2.res.gz
MEPFFDAYRAVKPFLVTSGNPPTKERIQSPTDRARYDDTTKCILCACCTTSCPVYWSEGSYFGPAAIVNAHRFIFDSRDEAAAERL--------------DILNEVDGVWRCRTTFNCTEACPRGIQVTQAIQEVKRALMFA


>d4f6ia_.pdb
--------------MV---NWAAVVDDFYQELFKAHPEYQNKFGFKGVALGSLKGNAAYKTQAGKTVDYINAAIGG-----SA--DAAGLASRHKGRNVGSAEFHNAKACLAKACSAHG-A---PD--LG-WAIDDILSHL----------
>d7diha_.pdb
MLSEETIRVIKSTVPLLKEHGTEITARMFELLFSKYPKTKELFAGA----------SE--EQPKKLANAIIAYATYIDRLEELDNAISTIARSHVRRNVKPEHYPLVKECLLQAIEEVLNPGEEVLKAWEEAYDFLAKTLITLEKKLYSQP




(base) D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output>D:\dummy\vscode_projects\cpp\TM-align_mass\bin\TMalign_mass.exe D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6iyia_.pdb D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d7diha_.pdb

 *********************************************************************
 * TM-align (Version 20210224): protein structure alignment          *
 * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005) *
 * Please email comments and suggestions to yangzhanglab@umich.edu   *
 *********************************************************************

Name of Chain_1: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6iyia_.pdb (to be superimposed onto Chain_2)
Name of Chain_2: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d7diha_.pdb
Length of Chain_1: 142 residues
Length of Chain_2: 139 residues

Aligned length= 129, RMSD=   1.93, Seq_ID=n_identical/n_aligned= 0.186
TM-score= 0.79753 (if normalized by length of Chain_1, i.e., LN=142, d0=4.43)
TM-score= 0.81284 (if normalized by length of Chain_2, i.e., LN=139, d0=4.38)
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)
SLTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSSA-QVHTHGKKIIDALTTAVNHID-DI-TG-ALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR------
:::::::::::::::::::::::::::::::::::::::::::::.         :  :::::::::::::::::: :. .: :::::::::: ::::::::::::::::::::::::: : ::::::::::::::::::::::::
MLSEETIRVIKSTVPLLKEHGTEITARMFELLFSKYPKTKELFAGA--------SE--EQPKKLANAIIAYATYIDRLEELDNAISTIARSHV-RRNVKPEHYPLVKECLLQAIEEVLN-P-GEEVLKAWEEAYDFLAKTLITLEKKLYSQP

Total CPU time is  0.02 seconds

(base) D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output>D:\dummy\vscode_projects\cpp\TM-align_mass\bin\TMalign_mass.exe D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6iyia_.pdb D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d4f6ia_.pdb

 *********************************************************************
 * TM-align (Version 20210224): protein structure alignment          *
 * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005) *
 * Please email comments and suggestions to yangzhanglab@umich.edu   *
 *********************************************************************

Name of Chain_1: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6iyia_.pdb (to be superimposed onto Chain_2)
Name of Chain_2: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d4f6ia_.pdb
Length of Chain_1: 142 residues
Length of Chain_2: 110 residues

Aligned length= 107, RMSD=   2.60, Seq_ID=n_identical/n_aligned= 0.112
TM-score= 0.59806 (if normalized by length of Chain_1, i.e., LN=142, d0=4.43)
TM-score= 0.73357 (if normalized by length of Chain_2, i.e., LN=110, d0=3.86)
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)
SLTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSS--A-QVHTHGKKIIDALTTAVNHIDDITGALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR
               ::  ::::::::::::::::::::::::..::::::::  : :::::::::::::::::::.  .  ::::::::::: :::::::::::::::::::::::   .   .:  ::::::::::: :
---------------MV--NWAAVVDDFYQELFKAHPEYQNKFGFKGVALGSLKGNAAYKTQAGKTVDYINAAIGGS--A--DAAGLASRHKG-RNVGSAEFHNAKACLAKACSAHG---A---PD--LGWAIDDILSH-L----

Total CPU time is  0.02 seconds

(base) D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output>D:\dummy\vscode_projects\cpp\TM-align_mass\bin\TMalign_mass.exe D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6iyia_.pdb D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6lumb2.pdb

 *********************************************************************
 * TM-align (Version 20210224): protein structure alignment          *
 * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005) *
 * Please email comments and suggestions to yangzhanglab@umich.edu   *
 *********************************************************************

Name of Chain_1: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6iyia_.pdb (to be superimposed onto Chain_2)
Name of Chain_2: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d6lumb2.pdb
Length of Chain_1: 142 residues
Length of Chain_2: 128 residues

Aligned length= 97, RMSD=   4.45, Seq_ID=n_identical/n_aligned= 0.103
TM-score= 0.40256 (if normalized by length of Chain_1, i.e., LN=142, d0=4.43)
TM-score= 0.43009 (if normalized by length of Chain_2, i.e., LN=128, d0=4.19)
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)
S------------------------------LTSADKSHVRSIWSKAGGSAEEIGAEALGRMLESFPNTKTYFDHYADLSVSSAQVHTHGKKIIDALTTAVNHIDDI-TGALSSLSTLHAQTLRVDPANFKILSHTILVVLALYFPADFTPEVHLACDKFLANVSHALADNYR
                               .....:: ::  :: :::: :::::::  ::..  .  :  :  .         .  ::  .:::::::::::::  ..::::::::.:  : :...:..:  :: ::::::::  ::   .: ..: ::::::::::..
-MEPFFDAYRAVKPFLVTSGNPPTKERIQSPTDRARYD-DT--TK-CILC-ACCTTSC--PVYW--S--E--G--S---------Y--FG--PAAIVNAHRFIFDS-RDEAAAERLDILN--E-VDGVWRCR--TT-FNCTEACP--RG---IQ-VTQ-AIQEVKRALMFA--

Total CPU time is  0.03 seconds

(base) D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output>


(base) D:\dummy\vscode_projects\matrix_align>D:\dummy\vscode_projects\cpp\TM-align_mass\bin\TMalign_mass.exe D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d4f6ia_.pdb D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d7diha_.pdb

 *********************************************************************
 * TM-align (Version 20210224): protein structure alignment          *
 * References: Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005) *
 * Please email comments and suggestions to yangzhanglab@umich.edu   *
 *********************************************************************

Name of Chain_1: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d4f6ia_.pdb (to be superimposed onto Chain_2)
Name of Chain_2: D:\dummy\vscode_projects\matrix_align\example_files\esm2_650m_example_output\d7diha_.pdb
Length of Chain_1: 110 residues
Length of Chain_2: 139 residues

Aligned length= 98, RMSD=   2.37, Seq_ID=n_identical/n_aligned= 0.173
TM-score= 0.72849 (if normalized by length of Chain_1, i.e., LN=110, d0=3.86)
TM-score= 0.59524 (if normalized by length of Chain_2, i.e., LN=139, d0=4.38)
(You should use TM-score normalized by length of the reference structure)

(":" denotes residue pairs of d <  5.0 Angstrom, "." denotes other aligned residues)
--------------MV---NWAAVVDDFYQELFKAHPEYQNKFGFKGVALGSLKGNAAYKTQAGKTVDYINAAIGG-----SA--DAAGLASRHKGRNVGSAEFHNAKACLAKACSAHG-A---PD--LG-WAIDDILSHL----------
              ::   :::::::::::::::::::::::::..          .:  ::::::::::::::::     ::  :::::::::::::::::::::::::::::::::: :   ::  :: ::::::::::
MLSEETIRVIKSTVPLLKEHGTEITARMFELLFSKYPKTKELFAGA----------SE--EQPKKLANAIIAYATYIDRLEELDNAISTIARSHVRRNVKPEHYPLVKECLLQAIEEVLNPGEEVLKAWEEAYDFLAKTLITLEKKLYSQP

Total CPU time is  0.00 seconds



esm2_650m_value_absdiff_swiss10000.dat
esm2_650m_value_absdiff_swiss10000.dat.stats

python scripts/esm_feature_check.py --infile uniprot_sprot.fasta.gz  --outfile esm_value_absdiff_swiss10000.dat  --random_seed 123 --device cuda --num_samples 10000 --model_path ＜ウエイトを保存したディレクトリへのパス＞/esm2_t33_650M_UR50D.pt --normalize True
として、断片化に影響を受ける値のインデクスをとった。
esm650m_value_absdiff_swiss10000.dat
の 2 カラム目の値が小さいほうが、配列の末端のアミノ酸が削れた際に値が変わりにくいもの。
cut -f 1 matrix_align/example_files/esm2_650m_example_output/esm2_650m_value_absdiff_swiss10000.dat |head -n 1024 > matrix_align/example_files/esm2_650m_example_output/esm2_650m_robust_index.dat
とし、トップ 1024 のインデクスをとった。

uniprot_sprot.fasta.gz は 20240523 に
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
UniProt: the universal protein knowledgebase in 2021. Nucleic acids research, 2021, 49.D1: D480-D489.
から DL した。




