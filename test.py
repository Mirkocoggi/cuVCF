import VCFParallelLib as vcf
import pandas as pd

res = vcf.vcf_parsed()
res.run("data/tiny.vcf", 1)
df1tmp = res.var_columns
df2tmp = res.alt_columns
df3tmp = res.samp_columns
df4tmp = res.alt_sample

#Abbiamo il problema della lunghezza degli array che deve essere uguale per tutti
data1 = {
    'var_number': df1tmp.var_number,
    'chrom': df1tmp.chrom,
    'pos': df1tmp.pos,
    'id':  df1tmp.id,
    'ref': df1tmp.ref,
    #'alt': df1tmp.alt,
    'qual': df1tmp.qual,
    'filter': df1tmp.filter
    #'info': df1tmp.info
}
#è necessario aggiungere dinamicamente le colonne deelle varie strutture
for elem in df1tmp.in_float:
    if len(elem.i_float) > 0:
        data1[elem.name] = elem.i_float
for elem in df1tmp.in_flag:
    if len(elem.i_flag) > 0:
        data1[elem.name] = elem.i_flag
for elem in df1tmp.in_string:
    if len(elem.i_string) > 0:
        data1[elem.name] = elem.i_string
for elem in df1tmp.in_int:
    if len(elem.i_int) > 0:
        data1[elem.name] = elem.i_int
df1 = pd.DataFrame(data1)

data2 = {
    'var_id': df2tmp.var_id,
    'alt_id': df2tmp.alt_id,
    'alt': df2tmp.alt
}
for elem in df2tmp.alt_float:
    if len(elem.i_float) > 0:
        data2[elem.name] = elem.i_float
for elem in df2tmp.alt_flag:
    if len(elem.i_flag) > 0:
        data2[elem.name] = elem.i_flag
for elem in df2tmp.alt_string:
    if len(elem.i_string) > 0:
        data2[elem.name] = elem.i_string
for elem in df2tmp.alt_int:
    if len(elem.i_int) > 0:
        data2[elem.name] = elem.i_int
df2 = pd.DataFrame(data2)

data3 = {
    'var_id': df3tmp.var_id,
    'samp_id': df3tmp.samp_id
}
for elem in df3tmp.samp_float:
    if len(elem.i_float) > 0:
        data3[elem.name] = elem.i_float
for elem in df3tmp.samp_flag:
    if len(elem.i_flag) > 0:
        data3[elem.name] = elem.i_flag
for elem in df3tmp.samp_string:
    if len(elem.i_string) > 0:
        data3[elem.name] = elem.i_string
for elem in df3tmp.samp_int:
    if len(elem.i_int) > 0:
        data3[elem.name] = elem.i_int
df3 = pd.DataFrame(data3)

data4 = {
    'var_id': df4tmp.var_id,
    'samp_id': df4tmp.samp_id,
    'alt_id': df4tmp.alt_id
}
for elem in df4tmp.samp_float:
    if len(elem.i_float) > 0:
        data4[elem.name] = elem.i_float
for elem in df4tmp.samp_flag:
    if len(elem.i_flag) > 0:
        data4[elem.name] = elem.i_flag
for elem in df4tmp.samp_string:
    if len(elem.i_string) > 0:
        data4[elem.name] = elem.i_string
for elem in df4tmp.samp_int:
    if len(elem.i_int) > 0:
        data4[elem.name] = elem.i_int
df4 = pd.DataFrame(data4)

print("DATAFRAME 1:")
print(df1)
print("DATAFRAME 2:")
print(df2)
print("DATAFRAME 3:")
print(df3)
print("DATAFRAME 4:")
print(df4)