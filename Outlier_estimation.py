# Importing necessary libraries for data manipulation, statistical analysis, and machine learning
import pandas as pd
import numpy as np
import random
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KernelDensity
from scipy import stats as st
import os

# Set the working directory path where the data files are located
working_path = ''

# Read and preprocess GWAS SNP data
print("Reading GWAS data")
gwas_filename = working_path + ""  # File path for the GWAS summary statistics
gwas_df = pd.read_csv(gwas_filename, sep='\t')  
# Convert 'pos' and 'chr' columns to integer type for consistency in data types
gwas_df['pos'] = gwas_df['pos'].astype(int)
gwas_df['chr'] = gwas_df['chr'].astype(int)
# Convert 'Beta' column to float for precise calculations
gwas_df['Beta'] = gwas_df['Beta'].astype(float)
# Extract and compute absolute values of the weights (Beta values)
weights = gwas_df['Beta'].to_numpy()
abs_weights = np.absolute(weights)

# Read and prepare control SNP data matched using SNPSnap
print("Reading controls")
cont_filename = working_path + ""  # File path for the SNPSnap matched control file
cont = pd.read_csv(cont_filename, sep='\t', header=None)  # Reading control data
# Fill missing values in the control data using forward fill method
cont.fillna(method='ffill', axis=0, inplace=True)
final_control_set = cont
# Remove the first column which contains "Input SNP" data, as it's not required for analysis
del final_control_set[final_control_set.columns[0]]

# Read iHS values from the 1000 Genomes Project (1KGP)
print("Reading iHS")
ihs_df = pd.read_csv('', sep='\t')  # File path for the iHS data file

# Merge GWAS SNP data with iHS values
print("Merge iHS and GWAS")
gwas_ihs_df = gwas_df.merge(ihs_df, on=['chr', 'pos'], how='left')  # Merging on chromosome and position
# Add absolute weights to the merged DataFrame for further analysis
gwas_ihs_df['abs_weights'] = abs_weights.tolist()
# Drop rows with missing values to ensure data integrity
gwas_ihs_df.dropna(inplace=True)

# List of populations to consider in the iHS analysis
ihs_pop_list = ["ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN", "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"]

# Setup for writing iHS AUC (Area Under Curve) results to a file
ihs_auc_file_name = working_path + ""  # Specify the file name for iHS AUC results
with open(ihs_auc_file_name, 'a+') as ot:
    # Write the list of populations as a header in the output file
    ot.write('\t'.join(ihs_pop_list) + '\n')

# Initialize bootstrap process for control sets
boot = 1000  # Number of bootstrap iterations to perform
bstat_auc = []  # List to store the bstat AUC values from each bootstrap iteration

# Looping through bootstrap iterations
for i in range(boot):
    print("Bootstrap:" + str(i))

    # Selecting a random column from the control set for the current bootstrap iteration
    rand1 = random.randrange(0, 3000, 1)  # Random index generation; adjust range based on the number of matched SNPs
    rand_arr = final_control_set.iloc[:, rand1].values  # Extracting a random set of control SNPs
    # Writing the randomly selected control SNPs to a file for record-keeping
    controls_file_name = working_path + ""  # File path for storing the selected control SNPs
    with open(controls_file_name, 'a+') as ot:
        ot.write('\t'.join(rand_arr) + '\n')

    # Parsing the SNP data into a list for further processing
    snp_list = []
    for j in range(len(rand_arr)):
        col = rand_arr[j].split(':')
        snp_list.append(col)

    # Preparing to calculate Bstat values for the control SNPs
    bstat_temp_file = working_path + "control_bed.txt"  # Temporary file to store control SNP data in BED format
    index_count = 0
    with open(bstat_temp_file, 'w') as f:
        for row in snp_list:
            f.write(row[0] + "\t" + row[1] + "\t" + str(index_count) + '\n')
            index_count += 1

    # Performing system calls to process the BED file and extract Bstat values
    f1 = working_path + "/temp/"
    # Convert to BED format, sort, and intersect with Bstat data using system commands
    com1 = "awk -v OFS='\t' '{print $1,$2-1,$2,$3}' " + f1 + "control_bed.txt > " + f1 + "control_bed.bed"
    os.system(com1)
    com2 = "sort " + f1 + "control_bed.bed > " + f1 + "control_bed_sorted.bed"
    os.system(com2)
    com3 = "bedtools intersect -a bstat_hg19.txt -b " + f1 + "control_bed_sorted.bed -wa -wb > " + f1 + "control_bstat.txt"  # Intersect with Bstat data
    os.system(com3)

    # Reading the processed Bstat values into a DataFrame
    snp_df = pd.read_csv('control_bstat.txt', sep='\t', header=None, index_col=None)
    # Merging Bstat data with GWAS data based on SNP indices
    snp_weight_match = pd.merge(snp_df, gwas_df, right_index=True, left_on=[7], how='left')

    # Extracting Bstat values and corresponding weights for KDE analysis
    x = snp_weight_match[3].to_numpy()
    wt = snp_weight_match['Beta'].to_numpy()
    wt_abs = np.absolute(wt)

    # Writing the Bstat values for each bootstrap iteration to a file
    bstat_file_name = working_path + "random_bootstrap_bstat_values.txt"
    with open(bstat_file_name, 'a+') as ot:
        ot.write('\t'.join(str(item) for item in x) + '\n')

    # Bstat AUC (Area Under the Curve) calculation
    print("Bstat")
    # Using GridSearchCV to find the optimal bandwidth for Kernel Density Estimation
    grid = GridSearchCV(KernelDensity(kernel='gaussian'), {'bandwidth': np.linspace(0.1, 2.0, 30)}, cv=5)
    grid.fit(x[:, None])
    best_bw = grid.best_params_
    # Calculating the standard deviation-adjusted bandwidth
    bw = best_bw['bandwidth'] / np.std(x)
    # Performing Kernel Density Estimation
    kde = st.gaussian_kde(x, bw_method=bw, weights=wt_abs)
    # Integrating the KDE to find the AUC
    ar = kde.integrate_box_1d(0, 317)
    # Appending the AUC result to a list
    bstat_auc.append(ar)

    # iHS analysis for controls
    print("iHS")
    # Preparing a DataFrame for iHS analysis
    snp_df = pd.DataFrame(snp_list, columns=['chr', 'pos'])
    snp_df['chr'] = snp_df['chr'].astype(int)
    snp_df['pos'] = snp_df['pos'].astype(float)

    # Merging the control SNP set with the iHS dataset
    merge_df = snp_df.merge(ihs_df, on=['chr', 'pos'], how='left')
    merge_df['abs_weights'] = abs_weights.tolist()
    merge_df.dropna(inplace=True)
    merge_df.drop(['chr', 'pos'], axis=1, inplace=True)
    # Extracting weights for the iHS analysis
    ihs_weights_control = merge_df['abs_weights'].to_numpy()

    # Initializing a list to store iHS AUC results
    ihs_auc = []
    for i in ihs_pop_list:  # Looping over all populations
        x = merge_df[i].to_numpy(dtype='float')
        # Skipping the population if it has less than 5 SNPs with iHS values
        if len(x) < 5:
            continue

        # Writing iHS values for each bootstrap iteration and population
        ihs_file_name = working_path + "/iHS_Bstat_Output/random_bootstrap_ihs_values_" + str(i) + ".txt"
        with open(ihs_file_name, 'a+') as ot:
            ot.write('\t'.join(str(item) for item in x) + '\n')

        # Estimating the bandwidth for Kernel Density Estimation
        grid = GridSearchCV(KernelDensity(kernel='gaussian'), {'bandwidth': np.linspace(0.1, 2.0, 30)}, cv=5)
        grid.fit(x[:, None])
        best_bw = grid.best_params_
        bw = best_bw['bandwidth'] / np.std(x)
        # Performing Kernel Density Estimation
        kde = st.gaussian_kde(x, bw_method=bw, weights=ihs_weights_control)
        # Calculating AUC for iHS
        ar = kde.integrate_box_1d(-2, 2)
        rev = 1 - ar
        ihs_auc.append(rev)

    # Writing the iHS AUC results for each population
    with open(ihs_auc_file_name, 'a+') as ot:
        ot.write('\t'.join(str(item) for item in ihs_auc) + '\n')

# Writing the calculated Bstat AUC values for each bootstrap iteration to a file
bstat_auc_file_name = working_path + "random_bootstrap_bstat_auc.txt"
with open(bstat_auc_file_name, 'a+') as ot:
    ot.write('\t'.join(str(item) for item in bstat_auc) + '\n')

# Processing GWAS Bstat and iHS data
f1 = working_path + "/temp/"
bstat_file_name = working_path + "random_bootstrap_bstat_values.txt"

# Extracting Bstat values for GWAS data
print("Get GWAS bstat")
# Preparing GWAS data in BED format for intersection with Bstat data
com1 = "awk -v OFS='\t' '{print $1,$2-1,$2,NR-1}' " + working_path + "/GWAS_clumped_hg19_bstat.txt > " + f1 + "gwas_bed.bed"
os.system(com1)
# Sorting the GWAS BED file
com2 = "sort " + f1 + "gwas_bed.bed > " + f1 + "gwas_bed_sorted.bed"
os.system(com2)
# Intersecting the sorted GWAS BED file with Bstat data to get corresponding Bstat values
com3 = "bedtools intersect -a bstat_hg19.txt -b " + f1 + "gwas_bed_sorted.bed -wa -wb > " + f1 + "gwas_bstat.txt"
os.system(com3)

# Reading the intersected Bstat values into a DataFrame
gwas_bstat_df = pd.read_csv('gwas_bstat.txt', sep='\t', header=None)
# Merging Bstat data with the original GWAS data
gwas_snp_weight_match = pd.merge(gwas_bstat_df, gwas_df, right_index=True, left_on=[7], how='left')
print(gwas_snp_weight_match)
# Dropping any rows with missing values
gwas_snp_weight_match.dropna(inplace=True)
# Extracting the Bstat values for further analysis
y = gwas_snp_weight_match[3].to_numpy(dtype='float')
print(y)
# Extracting the Beta weights for GWAS SNPs
gwas_wt = gwas_snp_weight_match['Beta'].to_numpy()
print(gwas_wt)
# Calculating the absolute values of the Beta weights
gwas_wt_abs = np.absolute(gwas_wt)

# Appending the extracted Bstat values to the bstat file
with open(bstat_file_name, 'a+') as ot:
    ot.write('\t'.join(str(item) for item in y) + '\n')

# Estimating the bandwidth for Kernel Density Estimation (KDE)
print("Get AUC")
grid = GridSearchCV(KernelDensity(kernel='gaussian'), {'bandwidth': np.linspace(0.1, 2.0, 30)}, cv=5)
# Fitting the grid search to the Bstat data
grid.fit(y[:, None])
# Extracting the best bandwidth parameter
best_bw = grid.best_params_
# Adjusting the bandwidth based on the standard deviation of the data
bw = best_bw['bandwidth'] / np.std(y)
# Performing KDE with the calculated bandwidth and weights
kde = st.gaussian_kde(y, bw_method=bw, weights=gwas_wt_abs)
# Calculating the Area Under the Curve (AUC) for the KDE
ar = kde.integrate_box_1d(0, 317)

# Reading previously calculated Bstat AUC values from a file
with open("random_bootstrap_bstat_auc.txt", "r") as file:
    for line in file:
        line = line.rstrip()
        bstat_au = line.split("\t")

# Converting the Bstat AUC values to a NumPy array for further processing
bstat_auc = np.array(bstat_au, dtype=float)


# Calculating the percentile of GWAS AUC (Area Under Curve) in comparison to the control Bstat AUC values
print("Get bstat percentile")
perc = st.percentileofscore(bstat_auc, ar, kind='strict')  # Computes the percentile score of 'ar' in the 'bstat_auc' array

# Writing the GWAS AUC and its percentile relative to controls to a file
bstat_result_file_name = working_path + "bstat_result.txt"  # File for storing the results
with open(bstat_result_file_name, 'a+') as ot:
    ot.write("GWAS AUC: " + str(ar) + '\n' + "Percentile wrt controls: " + str(perc) + "\n")

# Processing GWAS iHS data
print("Get GWAS iHS")
# Reading control AUC data for all populations
control_auc = pd.read_csv(ihs_auc_file_name, sep='\t')

# Looping through each population in the control AUC data
for column in control_auc:
    print("Get GWAS iHS AUC")
    # Dropping rows with missing values in the current population's column
    drop = gwas_ihs_df.dropna(subset=[column])
    # Extracting iHS values for the current population
    x = drop[column].to_numpy(dtype='float')
    print(x)
    # Calculating the absolute values of Beta weights for GWAS SNPs
    ihs_weights_gwas_abs = np.absolute(drop['Beta'].to_numpy())
    print(column)

    # Skipping populations with insufficient data
    if len(x) < 5:
        continue

    # Estimating the bandwidth for Kernel Density Estimation (KDE)
    grid = GridSearchCV(KernelDensity(kernel='gaussian'), {'bandwidth': np.linspace(0.1, 2.0, 30)}, cv=5)
    grid.fit(x[:, None])
    # Extracting the best bandwidth parameter
    best_bw = grid.best_params_
    bw = best_bw['bandwidth'] / np.std(x)
    # Performing KDE with the calculated bandwidth and weights
    kde = st.gaussian_kde(x, bw_method=bw, weights=ihs_weights_gwas_abs)
    # Calculating the AUC for iHS
    ar = kde.integrate_box_1d(-2, 2)
    rev = 1 - ar
    print(rev)

    # Extracting control AUC values for the current population
    c_auc = control_auc[column].to_numpy()
    # Calculating the percentile of GWAS AUC in comparison to control AUC values
    print("Get iHS percentile")
    perc = st.percentileofscore(c_auc, rev, kind='strict')
    print(perc)

    # Writing the GWAS iHS AUC and its percentile relative to controls for each population
    ihs_result_file_name = working_path + "ihs_result_" + str(column) + ".txt"
    with open(ihs_result_file_name, 'a+') as ot:
        ot.write("GWAS AUC: " + str(rev) + '\n' + "Percentile wrt controls: " + str(perc) + "\n")
