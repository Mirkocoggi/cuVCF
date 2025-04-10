import GPUParser as vcf
import pandas as pd
import time


res = vcf.vcf_parsed()

print("Start parsing")

start_run = time.perf_counter()
res.run("data/bos_taurus.vcf", 16)
end_run = time.perf_counter()
run_time = end_run - start_run
print(f"Tempo di esecuzione di run: {run_time:.4f} secondi")

start_df = time.perf_counter()
data1 = vcf.get_var_columns_data(res.var_columns)
data2 = vcf.get_alt_columns_data(res.alt_columns)
data3 = vcf.get_sample_columns_data(res.samp_columns)
data4 = vcf.get_alt_format_data(res.alt_sample)
end_df = time.perf_counter()
df_time = end_df - start_df
print(f"Tempo di creazione dei DataFrame c++: {df_time:.4f} secondi")

start_dfpd = time.perf_counter()
df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)
df3 = pd.DataFrame(data3)
df4 = pd.DataFrame(data4)
end_dfpd = time.perf_counter()
dfpd_time = end_dfpd - start_dfpd
print(f"Tempo di creazione dei DataFrame pandas: {dfpd_time:.4f} secondi")

start_csv = time.perf_counter()
df1.to_csv("df1.csv", index=False)
df2.to_csv("df2.csv", index=False)
df3.to_csv("df3.csv", index=False)
df4.to_csv("df4.csv", index=False)
end_csv = time.perf_counter()
csv_time = end_csv - start_csv
print(f"Tempo di salvataggio in CSV: {csv_time:.4f} secondi")

total_time = run_time + df_time + csv_time
print(f"Tempo totale: {total_time:.4f} secondi")

#print("DATAFRAME 1:")
#print(df1.iloc[:10])
#print("DATAFRAME 2:")
#print(df2.iloc[:10])
#print("DATAFRAME 3:")
#print(df3.iloc[:10])
#print("DATAFRAME 4:")
#print(df4.iloc[:10])

# Se necessario, in Python puoi reinterpretare gli array raw di half come np.float16:
#import numpy as np
# Ad esempio, per il campo "qual" in df1:
#if "qual" in data1:
#    # data1["qual"] è un array di uint16
#    arr = data1["qual"]
#    data1["qual"] = arr.view(np.float16)
#
#df1 = pd.DataFrame(data1)
#
#print("DATAFRAME 1.1:")
#print(df1.iloc[:10])