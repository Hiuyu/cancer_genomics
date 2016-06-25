function mutsig_perform(IN_MAF, OUT_PREFIX)
%IN_MAF='/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/2-mutsigcv_output/combine.maf'
%IN_MAF='/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/combine.annovar.maf';
%OUT_PREFIX='/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/2-mutsigcv_output/mutsig_output'
%OUT_PREFIX='/data/home/zuoxy/data/NPC/somatic/20160519_863closing/4-driver_gene/zz';
addpath('/data/home/zuoxy/apps/cancer_genomics/MutSigCV/MutSigCV_1.4');
MutSigCV(IN_MAF,'/data/home/zuoxy/apps/cancer_genomics/MutSigCV/dependencies/exome_full192.coverage.txt','/data/home/zuoxy/apps/cancer_genomics/MutSigCV/dependencies/gene.covariates.txt',OUT_PREFIX,'/data/home/zuoxy/apps/cancer_genomics/MutSigCV/dependencies/mutation_type_dictionary_file.txt','/data/home/zuoxy/apps/cancer_genomics/MutSigCV/dependencies/chr_files_hg19')
quit;
end
