import pandas as pd

input_filename_template = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/ase_aggregated_by_donor/{donor}.ase.lowthresh.{count_type}.phased.genelevel.tsv'

#subset donors

n_extra_donors = 3000

donor_list = ['HPSI0514i-puie_5','HPSI0214i-poih_4','HPSI0514i-letw_1','HPSI0813i-guss_1','HPSI0413i-nudd_1','HPSI1014i-sehl_6','HPSI0114i-joxm_1']
with open('/nfs/leia/research/stegle/dseaton/hipsci/singlecell_endodiff/data/list_of_singlecell_endodiff_donors.tsv', 'r') as f:
    extra_donors = [x.strip() for x in f.readlines()][:n_extra_donors]
    donor_list.extend(extra_donors)

donor_list = list(set(donor_list))[:]

output_filename_template = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/all_donors.ase.lowthresh.{count_type}.phased.genelevel.tsv'



#donor_list = ['HPSI1014i-sehl_6','HPSI0114i-joxm_1']
#output_filename_template = '/hps/nobackup/hipsci/scratch/singlecell_endodiff/data_processed/ase/test_subset_of_donors.ase.lowthresh.{count_type}.phased.genelevel.tsv'



# ase_df gives data for each gene, using phased SNP info to give proportion of expression from chrB

donor2cell_dict = dict()


for count_type in ['chrBcount','totalcount']:
    df_list = []
    for donor in donor_list[:]:
        try:
            filename = input_filename_template.format(donor=donor,count_type=count_type)
            df = pd.read_csv(filename,sep='\t',index_col=0)
#            df = df.reindex(gene_list)
            donor2cell_dict[donor] = list(df.columns)
            df_list.append(df)
        except:
            print('No file for cell line {}'.format(donor))
            pass
    print('Combining data from {} cell lines'.format(len(df_list)))
    df = pd.concat(df_list, axis=1)
    if count_type=='totalcount':
        total_df = df
        total_df.to_csv(output_filename_template.format(count_type=count_type), sep='\t')
    else:
        allelic_df = df
        allelic_df.to_csv(output_filename_template.format(count_type=count_type), sep='\t')

# donor list is reduced to whatever we actually have data for
donor_list = donor2cell_dict.keys()

