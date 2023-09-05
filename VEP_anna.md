**Данный** ~~текст~~ _отформатирован_ при помощи
# Создание VCF с дополнительными  альтернативными кодонами, кодирующими ту же аминокислоту. 

## 1. Подготовка аннотированного файла к п.3.


```python
import pandas as pd
import csv
import allel; print('scikit-allel', allel.__version__)
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
```

    scikit-allel 1.3.5
    


```python
file_name = "clinvar_vep_ref.vcf" # clinvar.vcf аннатированный VEP
clinvar_file_name = "clinvar.vcf" # https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
```

Создадим заголовок для датафрейма clinvar_ann_df


```python
def make_header(file_name):
    header_list = []
    with open(file_name, "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            if line.split()[0] == "#Uploaded_variation": # сохраним заголовок
                header_list.append(line.split())
                return header_list[0]
```


```python
header_list = make_header(file_name)
header_list
```




    ['#Uploaded_variation',
     'Location',
     'Allele',
     'Gene',
     'Feature',
     'Feature_type',
     'Consequence',
     'cDNA_position',
     'CDS_position',
     'Protein_position',
     'Amino_acids',
     'Codons',
     'Existing_variation',
     'Extra']



Загрузим clinvar.vcf


```python
%%time
df_clinvar = allel.vcf_to_dataframe(input = clinvar_file_name, fields='*')
df_clinvar
```

    CPU times: total: 15.6 s
    Wall time: 15.6 s
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CHROM</th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT_1</th>
      <th>ALT_2</th>
      <th>ALT_3</th>
      <th>QUAL</th>
      <th>AF_ESP</th>
      <th>AF_EXAC</th>
      <th>...</th>
      <th>MC</th>
      <th>ORIGIN</th>
      <th>RS</th>
      <th>SSR</th>
      <th>FILTER_PASS</th>
      <th>numalt</th>
      <th>altlen_1</th>
      <th>altlen_2</th>
      <th>altlen_3</th>
      <th>is_snp</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>925952</td>
      <td>1019397</td>
      <td>G</td>
      <td>A</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>SO:0001583|missense_variant</td>
      <td>1</td>
      <td>1640863258</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>925976</td>
      <td>1362713</td>
      <td>T</td>
      <td>C</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>SO:0001583|missense_variant</td>
      <td>1</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>926003</td>
      <td>1365270</td>
      <td>C</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>SO:0001583|missense_variant</td>
      <td>1</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>926014</td>
      <td>1377425</td>
      <td>G</td>
      <td>A</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>SO:0001575|splice_donor_variant</td>
      <td>1</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>930139</td>
      <td>1125147</td>
      <td>C</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>SO:0001627|intron_variant</td>
      <td>1</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>57893</th>
      <td>MT</td>
      <td>16159</td>
      <td>1525973</td>
      <td>C</td>
      <td>A</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>16</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>57894</th>
      <td>MT</td>
      <td>16179</td>
      <td>1525977</td>
      <td>CAA</td>
      <td>C</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>16</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>-2</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>57895</th>
      <td>MT</td>
      <td>16230</td>
      <td>1525975</td>
      <td>A</td>
      <td>G</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>16</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>57896</th>
      <td>MT</td>
      <td>16274</td>
      <td>1525974</td>
      <td>G</td>
      <td>A</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>16</td>
      <td>NaN</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>True</td>
    </tr>
    <tr>
      <th>57897</th>
      <td>NW_009646201.1</td>
      <td>83614</td>
      <td>17735</td>
      <td>TC</td>
      <td>T</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>...</td>
      <td>NaN</td>
      <td>1</td>
      <td>1556058284</td>
      <td>-1</td>
      <td>False</td>
      <td>1</td>
      <td>-1</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
<p>1303082 rows × 36 columns</p>
</div>



Создаем clinvar_ann_df из clinvar_vep_ref.vcf


```python
clinvar_ann_df = pd.read_table(file_name, header=None, skiprows=[i for i in range(0,49)])
clinvar_ann_df.columns = header_list
clinvar_ann_df
```

    C:\ProgramData\Miniconda3\lib\site-packages\IPython\core\interactiveshell.py:3369: DtypeWarning: Columns (3) have mixed types.Specify dtype option on import or set low_memory=False.
      exec(code_obj, self.user_global_ns, self.user_ns)
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>#Uploaded_variation</th>
      <th>Location</th>
      <th>Allele</th>
      <th>Gene</th>
      <th>Feature</th>
      <th>Feature_type</th>
      <th>Consequence</th>
      <th>cDNA_position</th>
      <th>CDS_position</th>
      <th>Protein_position</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>Existing_variation</th>
      <th>Extra</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1019397</td>
      <td>1:925952</td>
      <td>A</td>
      <td>148398</td>
      <td>NM_001385640.1</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>1057</td>
      <td>548</td>
      <td>183</td>
      <td>G/E</td>
      <td>gGg/gAg</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;GIVEN_REF=G;USED_REF=G</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1019397</td>
      <td>1:925952</td>
      <td>A</td>
      <td>148398</td>
      <td>NM_001385641.1</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>1057</td>
      <td>548</td>
      <td>183</td>
      <td>G/E</td>
      <td>gGg/gAg</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1019397</td>
      <td>1:925952</td>
      <td>A</td>
      <td>148398</td>
      <td>NM_152486.4</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>101</td>
      <td>11</td>
      <td>4</td>
      <td>G/E</td>
      <td>gGg/gAg</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;GIVEN_REF=G;USED_REF=G</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1019397</td>
      <td>1:925952</td>
      <td>A</td>
      <td>107985728</td>
      <td>NR_168405.1</td>
      <td>Transcript</td>
      <td>upstream_gene_variant</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>IMPACT=MODIFIER;DISTANCE=348;STRAND=-1;CANONIC...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1362713</td>
      <td>1:925976</td>
      <td>C</td>
      <td>148398</td>
      <td>NM_001385640.1</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>1081</td>
      <td>572</td>
      <td>191</td>
      <td>I/T</td>
      <td>aTc/aCc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;GIVEN_REF=T;USED_REF=T</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>12858667</th>
      <td>1525973</td>
      <td>MT:16159</td>
      <td>A</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>intergenic_variant</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>IMPACT=MODIFIER</td>
    </tr>
    <tr>
      <th>12858668</th>
      <td>1525977</td>
      <td>MT:16180-16181</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>intergenic_variant</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>IMPACT=MODIFIER</td>
    </tr>
    <tr>
      <th>12858669</th>
      <td>1525975</td>
      <td>MT:16230</td>
      <td>G</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>intergenic_variant</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>IMPACT=MODIFIER</td>
    </tr>
    <tr>
      <th>12858670</th>
      <td>1525974</td>
      <td>MT:16274</td>
      <td>A</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>intergenic_variant</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>IMPACT=MODIFIER</td>
    </tr>
    <tr>
      <th>12858671</th>
      <td>17735</td>
      <td>NW_009646201.1:83615</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>intergenic_variant</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>IMPACT=MODIFIER</td>
    </tr>
  </tbody>
</table>
<p>12858672 rows × 14 columns</p>
</div>



## 2. Вытаскиваем OMIM гены из аннотированного clinvar (clinvar_vep_ref.vcf)


```python
omim_list_file = 'omim_list.txt' # https://www.omim.org
```


```python
def omim_list_func(omim_list_file):
    omim_list = []
    with open(omim_list_file, "r") as file: #создаем список генов
        while True:
            line = file.readline()
            if not line:
                break
            omim_list.append(line)
    return omim_list
```

Удаляем NA и лишние символы


```python
omim_list = omim_list_func(omim_list_file)
omim_list_1 = [i.replace('\n', '') for i in omim_list ]
omim_list_2 = [i for i in omim_list_1 if 'NA' not in i]
omim_list_2
```




    ['NM_000014.6',
     'NM_144670.6',
     'NM_015665.6',
     'NM_024666.5',
     'NM_001605.3',
     'NM_020745.4',
     'NM_005763.4',
     'NM_020686.6',
     'NM_005502.4',
     'NM_173076.3',
     'NM_001089.3',
     'NM_000350.3',
     'NM_018672.5',
     'NM_019112.4',
     'NM_003742.4',
     'NM_000443.4',
     'NM_005689.4',
     'NM_004299.6',
     'NM_032583.4',
     'NM_000392.5',
     'NM_001171.6',
     'NM_001351295.2',
     'NM_000352.6',
     'NM_005691.4',
     'NM_000033.4',
     'NM_002858.4',
     'NM_005050.4',
     'NM_004827.3',
     'NM_022436.3',
     'NM_022437.3',
     'NM_001042472.3',
     'NM_016006.6',
     'NM_007313.2',
     'NM_198839.3',
     'NM_014049.5',
     'NM_000016.6',
     'NM_000017.4',
     'NM_001609.4',
     'NM_000018.4',
     'NM_013227.4',
     'NM_000019.4',
     'NM_001082486.2',
     'NM_001178057.2',
     'NM_000789.4',
     'NM_001300953.2',
     'NM_001122951.3',
     'NM_001098.3',
     'NM_004035.7',
     'NM_003500.4',
     'NM_001610.4',
     'NM_033068.3',
     'NM_001111035.3',
     'NM_004458.3',
     'NM_001318509.2',
     'NM_001100.4',
     'NM_001613.4',
     'NM_001101.5',
     'NM_005159.5',
     'NM_001614.5',
     'NM_001615.4',
     'NM_001130004.2',
     'NM_001103.4',
     'NM_001278344.2',
     'NM_004924.6',
     'NM_001105.5',
     'NM_000020.3',
     'NM_000666.3',
     'NM_000022.4',
     'NM_001282225.2',
     'NM_001110.4',
     'NM_003183.6',
     'NM_021723.5',
     'NM_030957.4',
     'NM_139025.5',
     'NM_139057.4',
     'NM_199355.4',
     'NM_014244.5',
     'NM_014243.3',
     'NM_014694.4',
     'NM_019032.6',
     'NM_001111.5',
     'NM_138422.4',
     'NM_021116.4',
     'NM_018417.6',
     'NM_004036.5',
     'NM_183357.3',
     'NM_015270.5',
     'NM_016824.5',
     'NM_013447.4',
     'NM_005682.7',
     'NM_198569.3',
     'NM_032119.4',
     'NM_000669.5',
     'NM_006721.4',
     'NM_001123.4',
     'NM_015339.5',
     'NM_000682.7',
     'NM_000024.6',
     'NM_000025.3',
     'NM_000026.4',
     'NM_199165.2',
     'NM_001129.5',
     'NM_002025.4',
     'NM_014423.4',
     'NM_006796.3',
     'NM_001134.3',
     'NM_000027.4',
     'NM_152336.4',
     'NM_021831.6',
     'NM_018238.4',
     'NM_000642.3',
     'NM_006412.4',
     'NM_003659.4',
     'NM_198576.4',
     'NM_001138.2',
     'NM_001384479.1',
     'NM_001330701.2',
     'NM_001286715.1',
     'NM_031850.3',
     'NM_000030.3',
     'NM_031900.4',
     'NM_000687.4',
     'NM_001029882.3',
     'NM_017651.5',
     'NM_001621.5',
     'NM_001622.4',
     'NM_020661.4',
     'NM_004208.4',
     'NM_004757.4',
     'NM_006303.4',
     'NM_003977.4',
     'NM_014336.5',
     'NM_000383.4',
     'NM_000476.3',
     'NM_013411.5',
     'NM_001625.4',
     'NM_152327.5',
     'NM_007202.4',
     'NM_005751.5',
     'NM_001354.6',
     'NM_001818.5',
     'NM_005989.4',
     'NM_005163.2',
     'NM_001626.6',
     'NM_005465.7',
     'NM_000031.6',
     'NM_000032.5',
     'NM_002860.4',
     'NM_000693.4',
     'NM_000690.4',
     'NM_000382.3',
     'NM_003748.4',
     'NM_001080.3',
     'NM_005589.4',
     'NM_001182.5',
     'NM_184041.5',
     'NM_000035.4',
     'NM_019109.5',
     'NM_001013620.4',
     'NM_001004127.3',
     'NM_024105.4',
     'NM_001099922.3',
     'NM_144988.4',
     'NM_033087.4',
     'NM_005787.6',
     'NM_013339.4',
     'NM_024079.5',
     'NM_024740.2',
     'NM_001378454.1',
     'NM_001139.3',
     'NM_000698.5',
     'NM_001165960.1',
     'NM_021628.3',
     'NM_020778.5',
     'NM_000478.6',
     'NM_020919.4',
     'NM_006492.3',
     'NM_021926.4',
     'NM_014324.6',
     'NM_016519.6',
     'NM_182680.1',
     'NM_152424.4',
     'NM_000479.5',
     'NM_020547.3',
     'NM_015365.3',
     'NM_030943.4',
     'NM_000036.3',
     'NM_001368809.2',
     'NM_001025389.2',
     'NM_000481.4',
     'NM_212557.4',
     'NM_014495.4',
     'NM_139314.3',
     'NM_000037.4',
     'NM_001142445.2',
     'NM_001148.6',
     'NM_001386175.1',
     'NM_020987.5',
     'NM_001149.4',
     'NM_054027.6',
     'NM_015114.3',
     'NM_001256182.2',
     'NM_013275.6',
     'NM_014915.3',
     'NM_173551.5',
     'NM_018685.5',
     'NM_018075.5',
     'NM_031418.4',
     'NM_213599.3',
     'NM_001025356.3',
     'NM_000216.4',
     'NM_032208.3',
     'NM_058172.6',
     'NM_007193.5',
     'NM_001154.4',
     'NM_001283.5',
     'NM_003916.5',
     'NM_001039569.2',
     'NM_004069.6',
     'NM_003664.5',
     'NM_004644.5',
     'NM_001261826.3',
     'NM_006594.5',
     'NM_007347.5',
     'NM_004722.4',
     'NM_007077.5',
     'NM_014855.3',
     'NM_000038.6',
     'NM_005883.3',
     'NM_153000.5',
     'NM_000039.3',
     'NM_001643.2',
     'NM_052968.5',
     'NM_000384.3',
     'NM_000483.5',
     'NM_000041.4',
     'NM_145660.2',
     'NM_001370595.2',
     'NM_000484.4',
     'NM_012096.3',
     'NM_000485.3',
     'NM_175073.3',
     'NM_000486.6',
     'NM_001651.4',
     'NM_001170.3',
     'NM_000044.6',
     'NM_001655.5',
     'NM_001024226.2',
     'NM_006420.3',
     'NM_000045.4',
     'NM_020754.4',
     'NM_001185077.3',
     'NM_014629.4',
     'NM_001130955.2',
     'NM_004723.4',
     'NM_015185.3',
     'NM_001173479.2',
     'NM_006015.6',
     'NM_001374820.1',
     'NM_152641.4',
     'NM_182896.3',
     'NM_012106.4',
     'NM_004311.4',
     'NM_177976.3',
     'NM_015161.3',
     'NM_018076.5',
     'NM_001288767.2',
     'NM_001105247.2',
     'NM_025139.6',
     'NM_014862.4',
     'NM_005720.4',
     'NM_000487.6',
     'NM_000046.5',
     'NM_198709.3',
     'NM_014960.5',
     'NM_022786.3',
     'NM_139058.3',
     'NM_001127505.3',
     'NM_177924.5',
     'NM_001198800.3',
     'NM_001198799.3',
     'NM_004316.4',
     'NM_018489.3',
     'NM_000048.4',
     'NM_133436.3',
     'NM_000049.4',
     'NM_004318.4',
     'NM_018136.5',
     'NM_000050.4',
     'NM_015338.6',
     'NM_018263.6',
     'NM_030632.3',
     'NM_032810.4',
     'NM_001170535.3',
     'NM_033064.5',
     'NM_007348.4',
     'NM_004849.4',
     'NM_004044.7',
     'NM_015915.5',
     'NM_015459.5',
     'NM_000051.4',
     'NM_001007026.2',
     'NM_145178.4',
     'NM_001010986.3',
     'NM_022089.4',
     'NM_001160233.2',
     'NM_000701.8',
     'NM_000702.4',
     'NM_152296.5',
     'NM_173201.5',
     'NM_001681.4',
     'NM_170665.4',
     'NM_001001331.4',
     'NM_001683.5',
     'NM_001001344.2',
     'NM_014382.5',
     'NM_001001937.2',
     'NM_001687.5',
     'NM_001183.6',
     'NM_005765.3',
     'NM_012463.4',
     'NM_020632.3',
     'NM_001690.4',
     'NM_001692.4',
     'NM_001693.4',
     'NM_001696.4',
     'NM_000052.7',
     'NM_000053.4',
     'NM_016529.6',
     'NM_005603.6',
     'NM_145691.4',
     'NM_001184.4',
     'NM_000489.6',
     'NM_000332.3',
     'NM_013236.4',
     'NM_002973.4',
     'NM_004993.6',
     'NM_000333.4',
     'NM_001698.3',
     'NM_198433.3',
     'NM_001015878.2',
     'NM_015570.4',
     'NM_000490.5',
     'NM_000054.6',
     'NM_004655.4',
     'NM_004048.4',
     'NM_152490.5',
     'NM_080605.4',
     'NM_012200.4',
     'NM_194318.4',
     'NM_001478.5',
     'NM_001497.4',
     'NM_007255.3',
     'NM_006876.3',
     'NM_015681.5',
     'NM_001243473.2',
     'NM_030578.4',
     'NM_001701.4',
     'NM_004281.4',
     'NM_001143985.1',
     'NM_004656.4',
     'NM_000465.4',
     'NM_001195306.2',
     'NM_024649.5',
     'NM_024685.4',
     'NM_152618.3',
     'NM_031885.5',
     'NM_033028.5',
     'NM_152384.3',
     'NM_176824.3',
     'NM_198428.3',
     'NM_005581.5',
     'NM_001139441.1',
     'NM_000709.4',
     'NM_183050.4',
     'NM_003921.5',
     'NM_022893.4',
     'NM_138576.4',
     'NM_017429.3',
     'NM_017745.6',
     'NM_004328.5',
     'NM_018429.3',
     'NM_001178020.3',
     'NM_004183.4',
     'NM_001195.5',
     'NM_003571.4',
     'NM_001711.6',
     'NM_001164405.2',
     'NM_030762.3',
     'NM_001080512.3',
     'NM_001003800.2',
     'NM_139343.3',
     'NM_004305.4',
     'NM_001715.3',
     'NM_000057.4',
     'NM_013314.4',
     'NM_212550.5',
     'NM_012388.4',
     'NM_000712.4',
     'NM_001199.4',
     'NM_006129.5',
     'NM_001200.4',
     'NM_001202.6',
     'NM_133468.5',
     'NM_004329.3',
     'NM_001203.3',
     'NM_001204.7',
     'NM_014753.4',
     'NM_212552.3',
     'NM_199186.3',
     'NM_004459.7',
     'NM_004333.6',
     'NM_152743.4',
     'NM_007294.4',
     'NM_000059.4',
     'NM_207189.4',
     'NM_001519.4',
     'NM_032043.3',
     'NM_001003694.2',
     'NM_153252.5',
     'NM_001130702.2',
     'NM_032667.6',
     'NM_001122955.4',
     'NM_057176.3',
     'NM_001370658.1',
     'NM_001287344.2',
     'NM_000061.3',
     'NM_001304561.1',
     'NM_001211.6',
     'NM_007073.4',
     'NM_020374.4',
     'NM_138425.4',
     'NM_152269.5',
     'NM_001130010.3',
     'NM_001031726.3',
     'NM_015991.4',
     'NM_000491.5',
     'NM_001212.4',
     'NM_172369.5',
     'NM_015645.5',
     'NM_001733.7',
     'NM_001354346.2',
     'NM_201442.4',
     'NM_000063.6',
     'NM_001286577.2',
     'NM_015531.6',
     'NM_000064.4',
     'NM_007293.3',
     'NM_001735.3',
     'NM_000562.3',
     'NM_000066.4',
     'NM_177965.4',
     'NM_018325.5',
     'NM_001218.5',
     'NM_000067.3',
     'NM_000717.5',
     'NM_001739.2',
     'NM_004056.6',
     'NM_016366.3',
     'NM_145200.5',
     'NM_000068.4',
     'NM_001127222.2',
     'NM_023035.3',
     'NM_001127221.2',
     'NM_000718.4',
     'NM_001129836.2',
     'NM_001129837.2',
     'NM_000719.7',
     'NM_001167625.2',
     'NM_001129829.2',
     'NM_001128839.3',
     'NM_000720.4',
     'NM_000721.4',
     'NM_005183.4',
     'NM_018896.5',
     'NM_198380.3',
     'NM_021098.3',
     'NM_000069.3',
     'NM_172364.5',
     'NM_000726.5',
     'NM_006078.5',
     'NM_004341.5',
     'NM_001742.4',
     'NM_001164737.2',
     'NM_006888.6',
     'NM_001743.6',
     'NM_015981.4',
     'NM_001220.5',
     'NM_015215.4',
     'NM_138793.4',
     'NM_001198868.2',
     'NM_000070.3',
     'NM_004055.5',
     'NM_032415.7',
     'NM_024110.4',
     'NM_052813.5',
     'NM_001013838.3',
     'NM_024537.4',
     'NM_004291.4',
     'NM_003688.3',
     'NM_001126055.2',
     'NM_032977.4',
     'NM_012114.3',
     'NM_001228.4',
     'NM_001080125.2',
     'NM_001231.5',
     'NM_001232.4',
     'NM_000388.4',
     'NM_001042440.5',
     'NM_053054.4',
     'NM_001172896.2',
     'NM_001753.5',
     'NM_033337.3',
     'NM_012232.6',
     'NM_005188.4',
     'NM_005142.3',
     'NM_000071.3',
     'NM_005189.3',
     'NM_017721.5',
     'NM_020785.2',
     'NM_001080522.2',
     'NM_133459.4',
     'NM_213607.3',
     'NM_144577.4',
     'NM_032357.4',
     'NM_016474.5',
     'NM_014008.5',
     'NM_024296.5',
     'NM_020198.3',
     'NM_178335.3',
     'NM_033124.5',
     'NM_001031737.3',
     'NM_032040.5',
     'NM_001135597.2',
     'NM_001080414.4',
     'NM_002986.3',
     'NM_002982.4',
     'NM_031443.4',
     'NM_003880.4',
     'NM_053056.3',
     'NM_001759.4',
     'NM_001099402.2',
     'NM_021147.5',
     'NM_152274.5',
     'NM_012073.5',
     'NM_006016.6',
     'NM_001770.6',
     'NM_198053.3',
     'NM_001242.5',
     'NM_001001547.3',
     'NM_000732.6',
     'NM_000733.4',
     'NM_000073.3',
     'NM_001250.6',
     'NM_000074.3',
     'NM_002389.4',
     'NM_000574.5',
     'NM_203330.2',
     'NM_001252.5',
     'NM_001783.4',
     'NM_000626.4',
     'NM_004356.4',
     'NM_001768.7',
     'NM_198196.3',
     'NM_138477.4',
     'NM_033312.3',
     'NM_001791.4',
     'NM_001178010.2',
     'NM_001254.4',
     'NM_024529.5',
     'NM_031942.5',
     'NM_004360.5',
     'NM_001797.4',
     'NM_004933.3',
     'NM_001171935.1',
     'NM_022124.6',
     'NM_001171933.1',
     'NM_001793.6',
     'NM_033100.4',
     'NM_052988.5',
     'NM_003718.5',
     'NM_000075.4',
     'NM_004935.4',
     'NM_018249.6',
     'NM_001259.8',
     'NM_001323289.2',
     'NM_001037343.2',
     'NM_003159.3',
     'NM_004064.5',
     'NM_000076.2',
     'NM_058195.4',
     'NM_000077.5',
     'NM_058197.5',
     'NM_016952.5',
     'NM_001264.5',
     'NM_030928.4',
     'NM_001039213.4',
     'NM_004364.5',
     'NM_001805.4',
     'NM_001807.6',
     'NM_001813.3',
     'NM_016343.4',
     'NM_018451.5',
     'NM_014704.4',
     'NM_153223.4',
     'NM_025009.5',
     'NM_014985.4',
     'NM_014956.5',
     'NM_032898.5',
     'NM_007186.6',
     'NM_025114.4',
     'NM_018718.3',
     'NM_001127182.2',
     'NM_014679.5',
     'NM_025180.5',
     'NM_001098802.3',
     'NM_016122.3',
     'NM_021267.5',
     'NM_178842.5',
     'NM_000078.3',
     'NM_021254.4',
     'NM_032930.3',
     'NM_004928.3',
     'NM_025145.7',
     'NM_001164496.2',
     'NM_145020.5',
     'NM_001039706.3',
     'NM_001710.6',
     'NM_032545.4',
     'NM_001928.4',
     'NM_000186.4',
     'NM_001014975.3',
     'NM_002113.3',
     'NM_021023.6',
     'NM_030787.4',
     'NM_000204.5',
     'NM_021914.8',
     'NM_002621.2',
     'NM_000492.4',
     'NM_001164144.3',
     'NM_001142933.2',
     'NM_020549.5',
     'NM_213720.3',
     'NM_001301339.2',
     'NM_016139.4',
     'NM_001270.4',
     'NM_001271.4',
     'NM_001005271.3',
     'NM_001273.5',
     'NM_017780.4',
     'NM_001316690.1',
     'NM_001170629.2',
     'NM_007194.4',
     'NM_001276.4',
     'NM_012110.4',
     'NM_003465.3',
     'NM_005198.5',
     'NM_000390.4',
     'NM_002768.5',
     'NM_001083314.4',
     'NM_014043.4',
     'NM_176812.5',
     'NM_001143981.2',
     'NM_000740.4',
     'NM_000079.4',
     'NM_001039523.3',
     'NM_000742.4',
     'NM_000744.7',
     'NM_000747.3',
     'NM_000751.3',
     'NM_000080.4',
     'NM_005199.5',
     'NM_018413.6',
     'NM_130468.4',
     'NM_004273.5',
     'NM_021615.5',
     'NM_001127896.2',
     'NM_014918.5',
     'NM_006384.4',
     'NM_006383.4',
     'NM_001304815.2',
     'NM_015125.5',
     'NM_022094.3',
     'NM_000246.4',
     'NM_001008388.5',
     'NM_001206999.2',
     'NM_006079.5',
     'NM_152515.5',
     'NM_013246.3',
     'NM_000083.3',
     'NM_004366.6',
     'NM_001830.4',
     'NM_000084.5',
     'NM_001287.6',
     'NM_004070.4',
     'NM_000085.5',
     'NM_021101.5',
     'NM_006984.5',
     'NM_144492.3',
     'NM_006580.4',
     'NM_148960.3',
     'NM_197947.3',
     'NM_001289.6',
     'NM_016929.5',
     'NM_001114086.2',
     'NM_024769.5',
     'NM_001042432.2',
     'NM_006493.4',
     'NM_017882.3',
     'NM_018941.4',
     'NM_006831.3',
     'NM_030813.6',
     'NM_006012.4',
     'NM_006660.5',
     'NM_174878.3',
     'NM_052995.2',
     'NM_004859.4',
     'NM_001127192.2',
     'NM_001298.3',
     'NM_001297.5',
     'NM_019098.5',
     'NM_017649.5',
     'NM_020184.4',
     'NM_001199303.2',
     'NM_006586.5',
     'NM_001843.4',
     'NM_005076.5',
     'NM_003632.3',
     'NM_014141.6',
     'NM_001008215.3',
     'NM_001012985.2',
     'NM_023077.3',
     'NM_025233.7',
     'NM_004086.3',
     'NM_018714.3',
     'NM_007357.3',
     'NM_015386.3',
     'NM_020751.3',
     'NM_153603.4',
     'NM_000493.4',
     'NM_001854.4',
     'NM_080680.3',
     'NM_004370.6',
     'NM_001130103.2',
     'NM_000494.4',
     'NM_001379500.1',
     'NM_000088.4',
     'NM_000089.4',
     'NM_198721.4',
     'NM_032888.4',
     'NM_001844.5',
     'NM_000090.4',
     'NM_001845.6',
     'NM_001846.4',
     'NM_000091.5',
     'NM_001130105.1',
     'NM_000092.5',
     'NM_000495.5',
     'NM_033380.3',
     'NM_001847.4',
     'NM_000093.5',
     'NM_000393.5',
     'NM_001848.3',
     'NM_058174.3',
     'NM_001849.4',
     'NM_004369.4',
     'NM_000094.4',
     'NM_005202.4',
     'NM_001851.6',
     'NM_001852.4',
     'NM_001853.4',
     'NM_006438.5',
     'NM_024027.5',
     'NM_024656.4',
     'NM_080538.2',
     'NM_005677.4',
     'NM_000095.3',
     'NM_000754.4',
     'NM_004371.4',
     'NM_004766.3',
     'NM_015697.9',
     'NM_016035.5',
     'NM_182476.3',
     'NM_182480.3',
     'NM_016138.5',
     'NM_020247.5',
     'NM_024876.4',
     'NM_020312.4',
     'NM_007074.4',
     'NM_001303.4',
     'NM_032901.4',
     'NM_001320976.2',
     'NM_004376.7',
     'NM_198076.6',
     'NM_032609.3',
     'NM_004373.4',
     'NM_001863.5',
     'NM_001866.3',
     'NM_004074.3',
     'NM_000096.4',
     'NM_020361.5',
     'NM_015692.5',
     'NM_023073.4',
     'NM_006651.4',
     'NM_001308.3',
     'NM_000097.7',
     'NM_001875.5',
     'NM_001876.4',
     'NM_001136052.3',
     'NM_000098.3',
     'NM_001006658.3',
     'NM_003805.5',
     'NM_000755.5',
     'NM_201253.3',
     'NM_173689.7',
     'NM_016302.4',
     'NM_052854.4',
     'NM_004380.3',
     'NM_001031717.4',
     'NM_015513.6',
     'NM_014171.6',
     'NM_004750.5',
     'NM_001101426.4',
     'NM_006371.5',
     'NM_000554.6',
     'NM_004075.5',
     'NM_000394.4',
     'NM_001885.3',
     'NM_005208.5',
     'NM_057093.2',
     'NM_001887.4',
     'NM_000496.3',
     'NM_004076.5',
     'NM_005210.4',
     'NM_020989.4',
     'NM_006891.4',
     'NM_017541.4',
     'NM_001888.5',
     'NM_005211.4',
     'NM_000395.3',
     'NM_172313.3',
     'NM_000760.4',
     'NM_001893.6',
     'NM_001895.4',
     'NM_024790.6',
     'NM_003476.5',
     'NM_000099.4',
     'NM_005213.4',
     'NM_000100.4',
     'NM_001328.3',
     'NM_025099.6',
     'NM_006565.4',
     'NM_004715.5',
     'NM_001902.6',
     'NM_005214.5',
     'NM_001903.5',
     'NM_001164883.2',
     'NM_004389.4',
     'NM_013266.4',
     'NM_001904.4',
     'NM_001085458.2',
     'NM_004937.3',
     'NM_001905.4',
     'NM_007272.3',
     'NM_000308.4',
     'NM_001908.5',
     'NM_001814.6',
     'NM_001909.5',
     'NM_003793.4',
     'NM_000396.4',
     'NM_001012759.3',
     'NM_001081.4',
     'NM_003590.5',
     'NM_003588.4',
     'NM_014780.5',
     'NM_001202543.2',
     'NM_181500.4',
     'NM_015267.4',
     'NM_005869.4',
     'NM_018294.6',
     'NM_001008540.2',
     'NM_003467.3',
     'NM_001915.4',
     'NM_001914.4',
     'NM_000398.7',
     'NM_000101.4',
     'NM_000397.4',
     'NM_001916.5',
     'NM_018947.6',
     'NM_001037333.3',
     'NM_015247.3',
     'NM_000497.4',
     'NM_000498.3',
     'NM_000102.4',
     'NM_031226.3',
     'NM_000104.4',
     'NM_000500.9',
     'NM_000782.5',
     'NM_183374.3',
     'NM_000784.4',
     'NM_000785.4',
     'NM_000762.6',
     'NM_000769.4',
     'NM_000771.4',
     'NM_000106.6',
     'NM_024514.5',
     'NM_183075.3',
     'NM_173483.4',
     'NM_207352.4',
     'NM_004820.5',
     'NM_152783.5',
     'NM_021080.5',
     'NM_016651.6',
     'NM_004393.6',
     'NM_172370.5',
     'NM_001349.4',
     'NM_018122.5',
     'NM_000787.4',
     'NM_001918.5',
     'NM_025000.4',
     'NM_015726.4',
     'NM_005215.4',
     'NM_016356.5',
     'NM_003737.4',
     'NM_001033855.3',
     'NM_001920.5',
     'NM_014026.6',
     'NM_004082.5',
     'NM_016286.4',
     'NM_000107.3',
     'NM_000790.4',
     'NM_001160147.2',
     'NM_001160148.2',
     'NM_015214.3',
     'NM_005216.5',
     'NM_006182.4',
     'NM_023935.3',
     'NM_152438.2',
     'NM_030653.4',
     'NM_001193416.3',
     'NM_016222.4',
     'NM_014314.4',
     'NM_001031725.6',
     'NM_021008.4',
     'NM_015213.4',
     'NM_001242896.3',
     'NM_001927.4',
     'NM_012079.6',
     'NM_003647.3',
     'NM_080916.3',
     'NM_014762.4',
     'NM_001360.3',
     'NM_024887.4',
     'NM_000791.4',
     'NM_021044.4',
     'NM_001361.5',
     'NM_018706.7',
     'NM_138615.3',
     'NM_014003.4',
     'NM_019887.6',
     'NM_005219.5',
     'NM_006729.5',
     'NM_001042517.2',
     'NM_177438.3',
     'NM_173602.3',
     'NM_152383.5',
     'NM_001363.5',
     'NM_001931.5',
     'NM_000108.5',
     'NM_020730.3',
     'NM_021120.4',
     'NM_016941.4',
     'NM_019074.4',
     'NM_005220.3',
     'NM_138281.3',
     'NM_005221.6',
     'NM_004006.3',
     'NM_000109.4',
     'NM_013391.3',
     'NM_004407.4',
     'NM_001081563.2',
     'NM_001174116.3',
     'NM_001080449.3',
     'NM_178452.6',
     'NM_001256714.1',
     'NM_130810.4',
     'NM_017802.4',
     'NM_015512.5',
     'NM_001277115.2',
     'NM_001372.4',
     'NM_012144.4',
     'NM_016306.6',
     'NM_153614.4',
     'NM_001039550.2',
     'NM_006736.6',
     'NM_058246.4',
     'NM_021800.3',
     'NM_145261.4',
     'NM_001012339.3',
     'NM_006260.5',
     ...]




```python
omim_pd = pd.DataFrame(omim_list_2, columns=['Feature'])
omim_pd
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Feature</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NM_000014.6</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NM_144670.6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NM_015665.6</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NM_024666.5</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NM_001605.3</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
    </tr>
    <tr>
      <th>3945</th>
      <td>NM_004773.4</td>
    </tr>
    <tr>
      <th>3946</th>
      <td>NM_207341.3</td>
    </tr>
    <tr>
      <th>3947</th>
      <td>NM_003460.2</td>
    </tr>
    <tr>
      <th>3948</th>
      <td>NM_001110354.2</td>
    </tr>
    <tr>
      <th>3949</th>
      <td>NM_020928.2</td>
    </tr>
  </tbody>
</table>
<p>3950 rows × 1 columns</p>
</div>



Найдем строчки в clinvar_vep_ref.vcf содержащие Feature из OMIM


```python
search_omim_df = clinvar_ann_df.merge(omim_pd, on='Feature', how='right')
search_omim_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>#Uploaded_variation</th>
      <th>Location</th>
      <th>Allele</th>
      <th>Gene</th>
      <th>Feature</th>
      <th>Feature_type</th>
      <th>Consequence</th>
      <th>cDNA_position</th>
      <th>CDS_position</th>
      <th>Protein_position</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>Existing_variation</th>
      <th>Extra</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>713107.0</td>
      <td>12:9074556</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>splice_donor_region_variant,intron_variant</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>-</td>
      <td>IMPACT=LOW;STRAND=-1;CANONICAL=YES;GIVEN_REF=C...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>781242.0</td>
      <td>12:9076832</td>
      <td>G</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>synonymous_variant</td>
      <td>3526</td>
      <td>3456</td>
      <td>1152</td>
      <td>Y</td>
      <td>taT/taC</td>
      <td>-</td>
      <td>IMPACT=LOW;STRAND=-1;CANONICAL=YES;GIVEN_REF=A...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>711716.0</td>
      <td>12:9079271</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3162</td>
      <td>3092</td>
      <td>1031</td>
      <td>R/Q</td>
      <td>cGa/cAa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>726982.0</td>
      <td>12:9079279</td>
      <td>G</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>synonymous_variant</td>
      <td>3154</td>
      <td>3084</td>
      <td>1028</td>
      <td>F</td>
      <td>ttT/ttC</td>
      <td>-</td>
      <td>IMPACT=LOW;STRAND=-1;CANONICAL=YES;GIVEN_REF=A...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>18171.0</td>
      <td>12:9079672</td>
      <td>C</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3068</td>
      <td>2998</td>
      <td>1000</td>
      <td>I/V</td>
      <td>Atc/Gtc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1503659</th>
      <td>155772.0</td>
      <td>5:61544156</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3502</td>
      <td>3487</td>
      <td>1163</td>
      <td>R/W</td>
      <td>Cgg/Tgg</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503660</th>
      <td>450504.0</td>
      <td>5:61544160</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3506</td>
      <td>3491</td>
      <td>1164</td>
      <td>H/R</td>
      <td>cAc/cGc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503661</th>
      <td>1424334.0</td>
      <td>5:61544165</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3511</td>
      <td>3496</td>
      <td>1166</td>
      <td>S/G</td>
      <td>Agt/Ggt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503662</th>
      <td>1417067.0</td>
      <td>5:61544201</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3547</td>
      <td>3532</td>
      <td>1178</td>
      <td>T/A</td>
      <td>Acc/Gcc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503663</th>
      <td>1353194.0</td>
      <td>5:61544238</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3584</td>
      <td>3569</td>
      <td>1190</td>
      <td>T/I</td>
      <td>aCa/aTa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
  </tbody>
</table>
<p>1503664 rows × 14 columns</p>
</div>



Выбираем только missense_variant и Canonical


```python
omim_missence_canonical_df = search_omim_df[(search_omim_df['Consequence'] == 'missense_variant') & (search_omim_df['Extra'].str.contains("CANONICAL=YES")==True)]
omim_missence_canonical_df['#Uploaded_variation']=omim_missence_canonical_df['#Uploaded_variation'].astype(int)
omim_missence_canonical_df
```

    C:\Users\Александр\AppData\Local\Temp\ipykernel_5964\3844229404.py:2: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      omim_missence_canonical_df['#Uploaded_variation']=omim_missence_canonical_df['#Uploaded_variation'].astype(int)
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>#Uploaded_variation</th>
      <th>Location</th>
      <th>Allele</th>
      <th>Gene</th>
      <th>Feature</th>
      <th>Feature_type</th>
      <th>Consequence</th>
      <th>cDNA_position</th>
      <th>CDS_position</th>
      <th>Protein_position</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>Existing_variation</th>
      <th>Extra</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2</th>
      <td>711716</td>
      <td>12:9079271</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3162</td>
      <td>3092</td>
      <td>1031</td>
      <td>R/Q</td>
      <td>cGa/cAa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>18171</td>
      <td>12:9079672</td>
      <td>C</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3068</td>
      <td>2998</td>
      <td>1000</td>
      <td>I/V</td>
      <td>Atc/Gtc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>5</th>
      <td>18172</td>
      <td>12:9079755</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>2985</td>
      <td>2915</td>
      <td>972</td>
      <td>C/Y</td>
      <td>tGt/tAt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>11</th>
      <td>18174</td>
      <td>12:9094987</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>2181</td>
      <td>2111</td>
      <td>704</td>
      <td>R/H</td>
      <td>cGt/cAt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>15</th>
      <td>729811</td>
      <td>12:9107574</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>899</td>
      <td>829</td>
      <td>277</td>
      <td>D/N</td>
      <td>Gac/Aac</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>1503659</th>
      <td>155772</td>
      <td>5:61544156</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3502</td>
      <td>3487</td>
      <td>1163</td>
      <td>R/W</td>
      <td>Cgg/Tgg</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503660</th>
      <td>450504</td>
      <td>5:61544160</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3506</td>
      <td>3491</td>
      <td>1164</td>
      <td>H/R</td>
      <td>cAc/cGc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503661</th>
      <td>1424334</td>
      <td>5:61544165</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3511</td>
      <td>3496</td>
      <td>1166</td>
      <td>S/G</td>
      <td>Agt/Ggt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503662</th>
      <td>1417067</td>
      <td>5:61544201</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3547</td>
      <td>3532</td>
      <td>1178</td>
      <td>T/A</td>
      <td>Acc/Gcc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>1503663</th>
      <td>1353194</td>
      <td>5:61544238</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3584</td>
      <td>3569</td>
      <td>1190</td>
      <td>T/I</td>
      <td>aCa/aTa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
  </tbody>
</table>
<p>432380 rows × 14 columns</p>
</div>




```python
df_1 = omim_missence_canonical_df.reset_index(drop=True)
df_1
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>#Uploaded_variation</th>
      <th>Location</th>
      <th>Allele</th>
      <th>Gene</th>
      <th>Feature</th>
      <th>Feature_type</th>
      <th>Consequence</th>
      <th>cDNA_position</th>
      <th>CDS_position</th>
      <th>Protein_position</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>Existing_variation</th>
      <th>Extra</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>711716</td>
      <td>12:9079271</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3162</td>
      <td>3092</td>
      <td>1031</td>
      <td>R/Q</td>
      <td>cGa/cAa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>18171</td>
      <td>12:9079672</td>
      <td>C</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3068</td>
      <td>2998</td>
      <td>1000</td>
      <td>I/V</td>
      <td>Atc/Gtc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>18172</td>
      <td>12:9079755</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>2985</td>
      <td>2915</td>
      <td>972</td>
      <td>C/Y</td>
      <td>tGt/tAt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>18174</td>
      <td>12:9094987</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>2181</td>
      <td>2111</td>
      <td>704</td>
      <td>R/H</td>
      <td>cGt/cAt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>729811</td>
      <td>12:9107574</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>899</td>
      <td>829</td>
      <td>277</td>
      <td>D/N</td>
      <td>Gac/Aac</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>432375</th>
      <td>155772</td>
      <td>5:61544156</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3502</td>
      <td>3487</td>
      <td>1163</td>
      <td>R/W</td>
      <td>Cgg/Tgg</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432376</th>
      <td>450504</td>
      <td>5:61544160</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3506</td>
      <td>3491</td>
      <td>1164</td>
      <td>H/R</td>
      <td>cAc/cGc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432377</th>
      <td>1424334</td>
      <td>5:61544165</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3511</td>
      <td>3496</td>
      <td>1166</td>
      <td>S/G</td>
      <td>Agt/Ggt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432378</th>
      <td>1417067</td>
      <td>5:61544201</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3547</td>
      <td>3532</td>
      <td>1178</td>
      <td>T/A</td>
      <td>Acc/Gcc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432379</th>
      <td>1353194</td>
      <td>5:61544238</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3584</td>
      <td>3569</td>
      <td>1190</td>
      <td>T/I</td>
      <td>aCa/aTa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
  </tbody>
</table>
<p>432380 rows × 14 columns</p>
</div>



Выберем нужные колонки #Uploaded_variation, Allele, Gene, cDNA_position, CDS_position, Protein_position, Amino_acids, Codons из df_1


```python
df_2 = df_1.iloc[:, [0,2,3,7,8,9,10,11]]
```

Создадим списки ref_amino_list, alt_amino_list, ref_codon_list, alt_codon_list


```python
def make_lists(df):
    for i, row in df.iterrows():
        ref_amino = row[6][:1]
        alt_amino = row[6][2:]
        ref_codon = row[7][:3]
        alt_codon = row[7][4:]
        yield ref_amino, alt_amino, ref_codon, alt_codon
```


```python
%%time
generator = make_lists(df_2)
all_list = []
ref_amino_list = []
alt_amino_list = []
ref_codon_list = []
alt_codon_list = []
i = 0
while True:
    try:
        all_list.append(next(generator))
        ref_amino_list.append(all_list[i][0])
        alt_amino_list.append(all_list[i][1])
        ref_codon_list.append(all_list[i][2])
        alt_codon_list.append(all_list[i][3])
        i+=1
    except StopIteration:
        break
```

    CPU times: total: 30.6 s
    Wall time: 30.6 s
    

## 3. Добавление альтернативного кодона, кодирующего ту же аминокислоту, отличающегося на один нуклеотид и которого не было в списке ранее

Создадим словари кодонов и аминокислот


```python
codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
}
```

Создадим функцию сравнения кодонов compare_triple. Если кодоны отличаются на одну букву, выводит референсный кодон, альтернативный кодон, референсный нуклеотид, альтернативный нуклеотид и позицию замены в кодоне


```python
def compare_triple(triple_a, triple_b):
    # Сравниваем два триплета, если различаются на одну букву, выводим triple_a, triple_b, a, b, position 
    list_dif_nucl = []
    pos = 0
    for x, y in zip(triple_a, triple_b):
        pos += 1
        if x != y:
            list_dif_nucl.append([x, y, pos])
    if len(list_dif_nucl) == 1:    
        return triple_a, triple_b, list_dif_nucl[0][0], list_dif_nucl[0][1], list_dif_nucl[0][2]
```

Создадим функцию сравнения референсного кодона с кортежем альтернативных кодонов ompare_ref_with_codons. Она посылает функции compare_triple, референсный кодон и альтернативный кодон из функции compare_lists. compare_ref_with_codons исключает альтернативного кодона, который уже был в clinvar_vep_ref.vcf. 


```python
def compare_ref_with_codons(ref_codon, alt_codons, alt_codon):
    list_compare_triple = []
    ref_pos = 1
    alt_pos = 1
    for i in alt_codons:
        list_a = []
        dif = 0
        if compare_triple(ref_codon, i) and i == alt_codon:
            ref_pos = compare_triple(ref_codon, i)[4]
    for i in alt_codons:
        list_a = []
        if compare_triple(ref_codon, i) and i != alt_codon:
            alt_pos = compare_triple(ref_codon, i)[4]
            list_a.append(compare_triple(ref_codon, i)[0:4])
            dif = alt_pos - ref_pos
            list_a.append(dif)    
            list_compare_triple.append(list_a)
    return list_compare_triple
```

Создадим функцию, которая создает список допольнительных альтернативных кодонов, кодирующих ту же аминокислоту, которых нет в исходном файле. На вход принимает списки референсных кодонов, альтернативной аминокислоты, альтернативных кодонов, которые получены из clinvar_vep_ref.vcf. Посылает функции compare_ref_with_codons референсный кодон, кортеж альтернативных кодонов, альтернативную аминокислотую.


```python
def compare_lists(ref_codon_list, alt_amino_list, alt_codon_list):
    additional_amino_list = []
    for ref_codon, alt_amino, alt_codon in zip(ref_codon_list, alt_amino_list, alt_codon_list):
        try:
            additional_amino_list.append(compare_ref_with_codons(ref_codon.upper(), codon_table[alt_amino], alt_codon.upper())) # сравниваем референсный кодон с кортежем альтернативных кодонов  
        except KeyError:
            additional_amino_list.append('#') # Если ошибка ключа помечаем "#"
            pass
    return additional_amino_list
```

Запускаем функцию compare_lists и создаем список additional_amino_list.


```python
additional_amino_list = compare_lists(ref_codon_list, alt_amino_list, alt_codon_list)
additional_amino_list
```




    [[],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAG', 'GAC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'CTG', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATA', 'G', 'A'), 0]],
     [],
     [],
     [],
     [],
     [[('AGT', 'CGT', 'A', 'C'), -2], [('AGT', 'AGA', 'T', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'TTG', 'A', 'T'), 0]],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [],
     [[('TGG', 'AGG', 'T', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [[('GGA', 'CGA', 'G', 'C'), 0]],
     [],
     [],
     [],
     [[('GGA', 'CGA', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [[('GGG', 'CGG', 'G', 'C'), 0]],
     [],
     [],
     [],
     [[('AGC', 'CGC', 'A', 'C'), -2], [('AGC', 'AGA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TTT', 'TTA', 'T', 'A'), 2], [('TTT', 'TTG', 'T', 'G'), 2]],
     [],
     [],
     [[('TGG', 'TGT', 'G', 'T'), 0]],
     [[('AAT', 'AAG', 'T', 'G'), 0]],
     [],
     [],
     [[('ATG', 'CTG', 'A', 'C'), 0]],
     [],
     [],
     [],
     [],
     [[('CAA', 'CAT', 'A', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAG', 'GAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TGG', 'AGG', 'T', 'A'), 0]],
     [],
     [],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [[('TTC', 'TTA', 'C', 'A'), 2], [('TTC', 'TTG', 'C', 'G'), 2]],
     [],
     [[('GTG', 'CTG', 'G', 'C'), 0]],
     [],
     '#',
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TTC', 'TTA', 'C', 'A'), 2], [('TTC', 'TTG', 'C', 'G'), 2]],
     [],
     [],
     [[('AAC', 'AAA', 'C', 'A'), 0]],
     [[('GGG', 'CGG', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGG', 'CGG', 'G', 'C'), 0]],
     [],
     [],
     [[('GAC', 'GAA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('TTT', 'TTA', 'T', 'A'), 2], [('TTT', 'TTG', 'T', 'G'), 2]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TGG', 'AGG', 'T', 'A'), 0]],
     [[('TGG', 'TGT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [],
     [],
     [[('AGT', 'AGA', 'T', 'A'), 2], [('AGT', 'AGG', 'T', 'G'), 2]],
     [],
     [],
     [],
     [[('CAA', 'CAC', 'A', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGG', 'AGG', 'G', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'TTG', 'A', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TTC', 'TTA', 'C', 'A'), 2], [('TTC', 'TTG', 'C', 'G'), 2]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TGG', 'AGG', 'T', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGG', 'CGG', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('AAC', 'AAA', 'C', 'A'), 0]],
     [[('CAG', 'CAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [[('GGA', 'CGA', 'G', 'C'), 0]],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [[('GAG', 'GAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGA', 'AGA', 'G', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('AGC', 'CGC', 'A', 'C'), -2], [('AGC', 'AGA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'CTG', 'G', 'C'), 0]],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGG', 'CGG', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAA', 'GAT', 'A', 'T'), 0]],
     [],
     [],
     [],
     [[('TGG', 'TGT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [[('TGC', 'AGC', 'T', 'A'), -1]],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TTC', 'TTG', 'C', 'G'), 0], [('TTC', 'CTC', 'T', 'C'), -2]],
     [],
     [[('TTT', 'TTA', 'T', 'A'), 0], [('TTT', 'CTT', 'T', 'C'), -2]],
     [[('GTG', 'CTG', 'G', 'C'), 0]],
     [[('GAT', 'GAG', 'T', 'G'), 0]],
     [],
     [],
     [],
     [[('CAC', 'CAA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('AAT', 'AAG', 'T', 'G'), 0]],
     [],
     [],
     [],
     [[('AAG', 'AAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'TTG', 'A', 'T'), 0]],
     [],
     [],
     [[('AAG', 'AAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('AAC', 'AAA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [],
     [[('GAT', 'GAG', 'T', 'G'), 0]],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [[('AAA', 'AAC', 'A', 'C'), 0]],
     [],
     [],
     [],
     [[('GGA', 'CGA', 'G', 'C'), 0]],
     [],
     [],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('TTC', 'TTG', 'C', 'G'), 0], [('TTC', 'CTC', 'T', 'C'), -2]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('AAG', 'AAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('AAC', 'AAA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [[('AGA', 'AGC', 'A', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('GAG', 'GAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAG', 'GAC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [[('TTC', 'TTG', 'C', 'G'), 0], [('TTC', 'CTC', 'T', 'C'), -2]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGG', 'CGG', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'CTG', 'G', 'C'), 0]],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATA', 'G', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGG', 'AGG', 'G', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAC', 'GAA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAT', 'G', 'T'), 0]],
     [],
     [],
     [[('GAC', 'GAA', 'C', 'A'), 0]],
     [[('GAC', 'GAA', 'C', 'A'), 0]],
     [[('CAG', 'CAC', 'G', 'C'), 0]],
     [[('AAC', 'AAA', 'C', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TTC', 'TTA', 'C', 'A'), 0], [('TTC', 'CTC', 'T', 'C'), -2]],
     [],
     [],
     [[('AAC', 'AAG', 'C', 'G'), 0]],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [[('GGG', 'CGG', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TTG', 'TTT', 'G', 'T'), 0]],
     [[('TTT', 'TTA', 'T', 'A'), 2], [('TTT', 'TTG', 'T', 'G'), 2]],
     [[('TGG', 'AGG', 'T', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAT', 'G', 'T'), 0]],
     [],
     [[('GGA', 'CGA', 'G', 'C'), 0]],
     [],
     [],
     [[('ATG', 'ATC', 'G', 'C'), 0], [('ATG', 'ATA', 'G', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'CTG', 'A', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'CTG', 'A', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TGG', 'AGG', 'T', 'A'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('AAT', 'AAG', 'T', 'G'), 0]],
     [],
     [],
     [],
     [],
     [[('CAG', 'CAC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GGG', 'AGG', 'G', 'A'), 0]],
     [[('AGG', 'AGC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [[('TTC', 'TTA', 'C', 'A'), 0], [('TTC', 'CTC', 'T', 'C'), -2]],
     [],
     [[('CAG', 'CAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [[('GTG', 'TTG', 'G', 'T'), 0]],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('TTG', 'TTT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATA', 'G', 'A'), 0]],
     [],
     [[('AGG', 'AGT', 'G', 'T'), 0]],
     [],
     [],
     [[('TGC', 'TCC', 'G', 'C'), 1]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ACC', 'TCC', 'A', 'T'), -1]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('GAG', 'GAT', 'G', 'T'), 0]],
     [[('TTT', 'TTA', 'T', 'A'), 2], [('TTT', 'TTG', 'T', 'G'), 2]],
     [[('GAG', 'GAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('AAG', 'AAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     [],
     [[('GGA', 'CGA', 'G', 'C'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('ATG', 'ATT', 'G', 'T'), 0], [('ATG', 'ATC', 'G', 'C'), 0]],
     [[('ACT', 'TCT', 'A', 'T'), -1]],
     [],
     [],
     [],
     [],
     [[('GAC', 'GAG', 'C', 'G'), 0]],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [],
     [[('AAG', 'AAT', 'G', 'T'), 0]],
     [],
     [],
     [],
     ...]



Из списка additional_amino_list создаем датафрейм additional_amino_df


```python
additional_amino_df = pd.DataFrame(additional_amino_list)
additional_amino_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>1</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>2</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>3</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>4</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>432375</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432376</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432377</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432378</th>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432379</th>
      <td>None</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
<p>432380 rows × 2 columns</p>
</div>



Соединяем колонки из df_2 (#Uploaded_variation, Allele, Gene, cDNA_position, CDS_position, Protein_position, Amino_acids, Codons; clinvar_vep_ref.vcf, missense_variant, Canonical) и additional_amino_df


```python
df_3 = pd.concat([df_2,additional_amino_df], axis=1)
df_3
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>#Uploaded_variation</th>
      <th>Allele</th>
      <th>Gene</th>
      <th>cDNA_position</th>
      <th>CDS_position</th>
      <th>Protein_position</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>0</th>
      <th>1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>711716</td>
      <td>T</td>
      <td>2</td>
      <td>3162</td>
      <td>3092</td>
      <td>1031</td>
      <td>R/Q</td>
      <td>cGa/cAa</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>1</th>
      <td>18171</td>
      <td>C</td>
      <td>2</td>
      <td>3068</td>
      <td>2998</td>
      <td>1000</td>
      <td>I/V</td>
      <td>Atc/Gtc</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>2</th>
      <td>18172</td>
      <td>T</td>
      <td>2</td>
      <td>2985</td>
      <td>2915</td>
      <td>972</td>
      <td>C/Y</td>
      <td>tGt/tAt</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>3</th>
      <td>18174</td>
      <td>T</td>
      <td>2</td>
      <td>2181</td>
      <td>2111</td>
      <td>704</td>
      <td>R/H</td>
      <td>cGt/cAt</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>4</th>
      <td>729811</td>
      <td>T</td>
      <td>2</td>
      <td>899</td>
      <td>829</td>
      <td>277</td>
      <td>D/N</td>
      <td>Gac/Aac</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>432375</th>
      <td>155772</td>
      <td>T</td>
      <td>57688</td>
      <td>3502</td>
      <td>3487</td>
      <td>1163</td>
      <td>R/W</td>
      <td>Cgg/Tgg</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432376</th>
      <td>450504</td>
      <td>G</td>
      <td>57688</td>
      <td>3506</td>
      <td>3491</td>
      <td>1164</td>
      <td>H/R</td>
      <td>cAc/cGc</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432377</th>
      <td>1424334</td>
      <td>G</td>
      <td>57688</td>
      <td>3511</td>
      <td>3496</td>
      <td>1166</td>
      <td>S/G</td>
      <td>Agt/Ggt</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432378</th>
      <td>1417067</td>
      <td>G</td>
      <td>57688</td>
      <td>3547</td>
      <td>3532</td>
      <td>1178</td>
      <td>T/A</td>
      <td>Acc/Gcc</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432379</th>
      <td>1353194</td>
      <td>T</td>
      <td>57688</td>
      <td>3584</td>
      <td>3569</td>
      <td>1190</td>
      <td>T/I</td>
      <td>aCa/aTa</td>
      <td>None</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
<p>432380 rows × 10 columns</p>
</div>



Создадим список нужных колонок Выберем нужные столбцы


```python
ncol = len(df_3.columns)
col_list = [0]
for i in range(6, ncol):
    col_list.append(i)
```

Выбираем нужные колонки


```python
df_4 = df_3.iloc[:, col_list]
df_4
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>#Uploaded_variation</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>0</th>
      <th>1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>711716</td>
      <td>R/Q</td>
      <td>cGa/cAa</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>1</th>
      <td>18171</td>
      <td>I/V</td>
      <td>Atc/Gtc</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>2</th>
      <td>18172</td>
      <td>C/Y</td>
      <td>tGt/tAt</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>3</th>
      <td>18174</td>
      <td>R/H</td>
      <td>cGt/cAt</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>4</th>
      <td>729811</td>
      <td>D/N</td>
      <td>Gac/Aac</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>432375</th>
      <td>155772</td>
      <td>R/W</td>
      <td>Cgg/Tgg</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432376</th>
      <td>450504</td>
      <td>H/R</td>
      <td>cAc/cGc</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432377</th>
      <td>1424334</td>
      <td>S/G</td>
      <td>Agt/Ggt</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432378</th>
      <td>1417067</td>
      <td>T/A</td>
      <td>Acc/Gcc</td>
      <td>None</td>
      <td>None</td>
    </tr>
    <tr>
      <th>432379</th>
      <td>1353194</td>
      <td>T/I</td>
      <td>aCa/aTa</td>
      <td>None</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
<p>432380 rows × 5 columns</p>
</div>



Выберем только строки, содержащие новые альтернативные кодоны


```python
df_5 = df_4.loc[df_4[0].isna()!=True]
df_5.rename(columns={'#Uploaded_variation': 'ID'}, inplace=True)
df_5.head(20)
```

    C:\ProgramData\Miniconda3\lib\site-packages\pandas\core\frame.py:5039: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      return super().rename(
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>ID</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>0</th>
      <th>1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>24</th>
      <td>372642</td>
      <td>E/D</td>
      <td>gaG/gaT</td>
      <td>[(GAG, GAC, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>30</th>
      <td>800252</td>
      <td>V/L</td>
      <td>Gtg/Ttg</td>
      <td>[(GTG, CTG, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>38</th>
      <td>644546</td>
      <td>Q/H</td>
      <td>caG/caT</td>
      <td>[(CAG, CAC, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>51</th>
      <td>1515743</td>
      <td>M/I</td>
      <td>atG/atC</td>
      <td>[(ATG, ATT, G, T), 0]</td>
      <td>[(ATG, ATA, G, A), 0]</td>
    </tr>
    <tr>
      <th>56</th>
      <td>561776</td>
      <td>S/R</td>
      <td>agT/agG</td>
      <td>[(AGT, CGT, A, C), -2]</td>
      <td>[(AGT, AGA, T, A), 0]</td>
    </tr>
    <tr>
      <th>65</th>
      <td>1438772</td>
      <td>M/L</td>
      <td>Atg/Ctg</td>
      <td>[(ATG, TTG, A, T), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>67</th>
      <td>1043900</td>
      <td>V/L</td>
      <td>Gtg/Ctg</td>
      <td>[(GTG, TTG, G, T), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>70</th>
      <td>860781</td>
      <td>W/R</td>
      <td>Tgg/Cgg</td>
      <td>[(TGG, AGG, T, A), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>76</th>
      <td>120256</td>
      <td>D/E</td>
      <td>gaC/gaA</td>
      <td>[(GAC, GAG, C, G), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>78</th>
      <td>1403638</td>
      <td>G/R</td>
      <td>Gga/Aga</td>
      <td>[(GGA, CGA, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>82</th>
      <td>942534</td>
      <td>G/R</td>
      <td>Gga/Aga</td>
      <td>[(GGA, CGA, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>87</th>
      <td>963431</td>
      <td>G/R</td>
      <td>Ggg/Agg</td>
      <td>[(GGG, CGG, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>91</th>
      <td>1019847</td>
      <td>S/R</td>
      <td>agC/agG</td>
      <td>[(AGC, CGC, A, C), -2]</td>
      <td>[(AGC, AGA, C, A), 0]</td>
    </tr>
    <tr>
      <th>100</th>
      <td>546698</td>
      <td>F/L</td>
      <td>Ttt/Ctt</td>
      <td>[(TTT, TTA, T, A), 2]</td>
      <td>[(TTT, TTG, T, G), 2]</td>
    </tr>
    <tr>
      <th>103</th>
      <td>939899</td>
      <td>W/C</td>
      <td>tgG/tgC</td>
      <td>[(TGG, TGT, G, T), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>104</th>
      <td>1008086</td>
      <td>N/K</td>
      <td>aaT/aaA</td>
      <td>[(AAT, AAG, T, G), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>107</th>
      <td>451644</td>
      <td>M/L</td>
      <td>Atg/Ttg</td>
      <td>[(ATG, CTG, A, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>112</th>
      <td>617557</td>
      <td>Q/H</td>
      <td>caA/caC</td>
      <td>[(CAA, CAT, A, T), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>146</th>
      <td>1465244</td>
      <td>E/D</td>
      <td>gaG/gaC</td>
      <td>[(GAG, GAT, G, T), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>150</th>
      <td>1403325</td>
      <td>M/I</td>
      <td>atG/atA</td>
      <td>[(ATG, ATT, G, T), 0]</td>
      <td>[(ATG, ATC, G, C), 0]</td>
    </tr>
  </tbody>
</table>
</div>



## 4. Создание VCF 

Свяжем df_5 и df_clinvar_2 (CHROM, POS, ID, REF, ALT_1, CLNSIG из clinvar.vcf)


```python
df_clinvar_2 = df_clinvar.iloc[:, [0,1,2,3,4,18]]
df_clinvar_2['ID']=df_clinvar_2['ID'].astype(int)
df_add = df_clinvar_2.merge(df_5, on='ID', how='right')
df_add
```

    C:\Users\Александр\AppData\Local\Temp\ipykernel_5964\2067860276.py:2: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      df_clinvar_2['ID']=df_clinvar_2['ID'].astype(int)
    




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CHROM</th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT_1</th>
      <th>CLNSIG</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>0</th>
      <th>1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>12</td>
      <td>8823290</td>
      <td>372642</td>
      <td>G</td>
      <td>T</td>
      <td>Uncertain_significance</td>
      <td>E/D</td>
      <td>gaG/gaT</td>
      <td>[(GAG, GAC, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>1</th>
      <td>12</td>
      <td>8823750</td>
      <td>800252</td>
      <td>G</td>
      <td>T</td>
      <td>Likely_benign</td>
      <td>V/L</td>
      <td>Gtg/Ttg</td>
      <td>[(GTG, CTG, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>2</th>
      <td>12</td>
      <td>8823848</td>
      <td>644546</td>
      <td>G</td>
      <td>T</td>
      <td>Uncertain_significance</td>
      <td>Q/H</td>
      <td>caG/caT</td>
      <td>[(CAG, CAC, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>3</th>
      <td>12</td>
      <td>8835599</td>
      <td>1515743</td>
      <td>G</td>
      <td>C</td>
      <td>Uncertain_significance</td>
      <td>M/I</td>
      <td>atG/atC</td>
      <td>[(ATG, ATT, G, T), 0]</td>
      <td>[(ATG, ATA, G, A), 0]</td>
    </tr>
    <tr>
      <th>4</th>
      <td>12</td>
      <td>8835653</td>
      <td>561776</td>
      <td>T</td>
      <td>G</td>
      <td>Uncertain_significance</td>
      <td>S/R</td>
      <td>agT/agG</td>
      <td>[(AGT, CGT, A, C), -2]</td>
      <td>[(AGT, AGA, T, A), 0]</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>53846</th>
      <td>5</td>
      <td>61521384</td>
      <td>1381845</td>
      <td>T</td>
      <td>G</td>
      <td>Uncertain_significance</td>
      <td>H/Q</td>
      <td>caT/caG</td>
      <td>[(CAT, CAA, T, A), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>53847</th>
      <td>5</td>
      <td>61538904</td>
      <td>1311009</td>
      <td>G</td>
      <td>C</td>
      <td>Uncertain_significance</td>
      <td>W/C</td>
      <td>tgG/tgC</td>
      <td>[(TGG, TGT, G, T), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>53848</th>
      <td>5</td>
      <td>61539627</td>
      <td>1347650</td>
      <td>A</td>
      <td>C</td>
      <td>Uncertain_significance</td>
      <td>E/D</td>
      <td>gaA/gaC</td>
      <td>[(GAA, GAT, A, T), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>53849</th>
      <td>5</td>
      <td>61543813</td>
      <td>1514466</td>
      <td>G</td>
      <td>T</td>
      <td>Uncertain_significance</td>
      <td>K/N</td>
      <td>aaG/aaT</td>
      <td>[(AAG, AAC, G, C), 0]</td>
      <td>None</td>
    </tr>
    <tr>
      <th>53850</th>
      <td>5</td>
      <td>61543957</td>
      <td>1490854</td>
      <td>G</td>
      <td>C</td>
      <td>Uncertain_significance</td>
      <td>K/N</td>
      <td>aaG/aaC</td>
      <td>[(AAG, AAT, G, T), 0]</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
<p>53851 rows × 10 columns</p>
</div>



Формируем VCF формат. Для этого создадим функцию, которая создает предварительный список, для датафрейма. Этот список содержит колонки для формата VCF: #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO


```python
def make_df_list(df):
    df_list=[]
    for row in df.itertuples():        
        for i in range(9, len(row)): 
            if row[i] and row[i] != '#':
                str_list = []
                pos = int(row[2]) + int(row[i][1])
                ID = str(row[3]) + '_' + str(i-8)
                REF = row[i][0][2]
                ALT_1 = row[i][0][3]
                REF_codon = row[i][0][0]
                ALT_codon = row[i][0][1]
                str_list.append(row[1])
                str_list.append(pos)
                str_list.append(ID)
                str_list.append(REF)
                str_list.append(ALT_1)
                str_list.append('.') # 'QUAL'
                str_list.append('.') # 'FILTER'
                info_str = ''
                for k in range(6, 8):
                    info_str += str(row[k]) + ';'
                info_str += REF_codon + ';' + ALT_codon
                str_list.append(info_str)
                df_list.append(str_list)
    return df_list
```


```python
df_list = make_df_list(df_add)
```

Создаем из списка df_list датафрейм df_6


```python
df_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
df_6 = pd.DataFrame(df_list, columns=df_columns)
df_6.set_index('#CHROM', inplace=True)
df_6
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT</th>
      <th>QUAL</th>
      <th>FILTER</th>
      <th>INFO</th>
    </tr>
    <tr>
      <th>#CHROM</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>12</th>
      <td>8823290</td>
      <td>372642_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;E/D;GAG;GAC</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8823750</td>
      <td>800252_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;V/L;GTG;CTG</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8823848</td>
      <td>644546_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;Q/H;CAG;CAC</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8835599</td>
      <td>1515743_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;M/I;ATG;ATT</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8835599</td>
      <td>1515743_2</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;M/I;ATG;ATA</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>5</th>
      <td>61521384</td>
      <td>1381845_1</td>
      <td>T</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;H/Q;CAT;CAA</td>
    </tr>
    <tr>
      <th>5</th>
      <td>61538904</td>
      <td>1311009_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;W/C;TGG;TGT</td>
    </tr>
    <tr>
      <th>5</th>
      <td>61539627</td>
      <td>1347650_1</td>
      <td>A</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;E/D;GAA;GAT</td>
    </tr>
    <tr>
      <th>5</th>
      <td>61543813</td>
      <td>1514466_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;K/N;AAG;AAC</td>
    </tr>
    <tr>
      <th>5</th>
      <td>61543957</td>
      <td>1490854_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Uncertain_significance;K/N;AAG;AAT</td>
    </tr>
  </tbody>
</table>
<p>63383 rows × 7 columns</p>
</div>



Сохраним df_6 в vcf. Для этого создадим функцию save_VCF


```python
def save_VCF(df, name):
    #надо исправить
    header = """##fileformat=
##fileDate=
##source=
##reference=
#CHROM POS ID REF ALT QUAL FILTER INFO
"""
    with open(name, 'w') as vcf:
        vcf.write(header)
    df.to_csv(name, sep="\t", mode='a', header=False)
```


```python
save_VCF(df_6, "Omim_long.vcf")
```

## 5. Создаем 4 датафрейма

Создадим 4 датафрейма benign_df, likely_benign_df, pathogenic_df, likely_pathogenic_df из df_26


```python
benign_df = df_6[(df_6['INFO'].str.contains("Benign")==True) & (df_6['INFO'].str.contains("Likely_benign")==False)]
benign_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT</th>
      <th>QUAL</th>
      <th>FILTER</th>
      <th>INFO</th>
    </tr>
    <tr>
      <th>#CHROM</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>12</th>
      <td>8838341</td>
      <td>120256_1</td>
      <td>C</td>
      <td>G</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;D/E;GAC;GAG</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8849672</td>
      <td>561635_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;V/L;GTG;CTG</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8852296</td>
      <td>1169145_1</td>
      <td>C</td>
      <td>G</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;D/E;GAC;GAG</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8867913</td>
      <td>413827_1</td>
      <td>C</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;N/K;AAC;AAA</td>
    </tr>
    <tr>
      <th>6</th>
      <td>44310332</td>
      <td>136225_1</td>
      <td>C</td>
      <td>G</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;D/E;GAC;GAG</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>16</th>
      <td>88434626</td>
      <td>320967_1</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;G/R;GGG;AGG</td>
    </tr>
    <tr>
      <th>1</th>
      <td>151286637</td>
      <td>779222_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;G/R;GGG;CGG</td>
    </tr>
    <tr>
      <th>1</th>
      <td>151289257</td>
      <td>717811_1</td>
      <td>A</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;Q/H;CAA;CAT</td>
    </tr>
    <tr>
      <th>7</th>
      <td>76425072</td>
      <td>791452_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;E/D;GAG;GAT</td>
    </tr>
    <tr>
      <th>7</th>
      <td>76441894</td>
      <td>712658_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Benign;R/S;AGG;AGC</td>
    </tr>
  </tbody>
</table>
<p>1510 rows × 7 columns</p>
</div>




```python
likely_benign_df = df_6[df_6['INFO'].str.contains("Likely_benign")==True]
likely_benign_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT</th>
      <th>QUAL</th>
      <th>FILTER</th>
      <th>INFO</th>
    </tr>
    <tr>
      <th>#CHROM</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>12</th>
      <td>8823750</td>
      <td>800252_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;V/L;GTG;CTG</td>
    </tr>
    <tr>
      <th>12</th>
      <td>8868027</td>
      <td>241909_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Benign/Likely_benign;E/D;GAG;GAT</td>
    </tr>
    <tr>
      <th>12</th>
      <td>53307533</td>
      <td>309718_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Benign/Likely_benign;G/R;GGG;CGG</td>
    </tr>
    <tr>
      <th>6</th>
      <td>44300553</td>
      <td>1201845_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;Q/H;CAG;CAT</td>
    </tr>
    <tr>
      <th>6</th>
      <td>44306988</td>
      <td>357073_1</td>
      <td>A</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;M/L;ATG;CTG</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>16</th>
      <td>88436510</td>
      <td>808133_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;G/R;GGA;CGA</td>
    </tr>
    <tr>
      <th>16</th>
      <td>88437796</td>
      <td>511466_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;R/S;AGG;AGC</td>
    </tr>
    <tr>
      <th>16</th>
      <td>88437797</td>
      <td>388236_1</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;G/R;GGG;AGG</td>
    </tr>
    <tr>
      <th>16</th>
      <td>88438445</td>
      <td>195113_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Benign/Likely_benign;G/R;GGG;CGG</td>
    </tr>
    <tr>
      <th>16</th>
      <td>88438876</td>
      <td>392227_1</td>
      <td>C</td>
      <td>G</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_benign;H/Q;CAC;CAG</td>
    </tr>
  </tbody>
</table>
<p>2367 rows × 7 columns</p>
</div>




```python
pathogenic_df = df_6[(df_6['INFO'].str.contains("Pathogenic")==True) & (df_6['INFO'].str.contains("Likely_pathogenic")==False)]
pathogenic_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT</th>
      <th>QUAL</th>
      <th>FILTER</th>
      <th>INFO</th>
    </tr>
    <tr>
      <th>#CHROM</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2</th>
      <td>214990791</td>
      <td>2863_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;G/R;GGA;CGA</td>
    </tr>
    <tr>
      <th>1</th>
      <td>94007725</td>
      <td>99418_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;G/R;GGA;CGA</td>
    </tr>
    <tr>
      <th>1</th>
      <td>94010868</td>
      <td>1184497_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;M/I;ATG;ATC</td>
    </tr>
    <tr>
      <th>1</th>
      <td>94010868</td>
      <td>1184497_2</td>
      <td>G</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;M/I;ATG;ATA</td>
    </tr>
    <tr>
      <th>1</th>
      <td>94019708</td>
      <td>1456885_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;M/I;ATG;ATT</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>8</th>
      <td>105801286</td>
      <td>156583_1</td>
      <td>A</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;S/R;AGT;CGT</td>
    </tr>
    <tr>
      <th>8</th>
      <td>105801288</td>
      <td>156583_2</td>
      <td>T</td>
      <td>G</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;S/R;AGT;AGG</td>
    </tr>
    <tr>
      <th>X</th>
      <td>137567449</td>
      <td>545553_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;C/S;TGC;TCC</td>
    </tr>
    <tr>
      <th>X</th>
      <td>137567448</td>
      <td>11437_1</td>
      <td>T</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;C/S;TGC;AGC</td>
    </tr>
    <tr>
      <th>1</th>
      <td>40285988</td>
      <td>4272_1</td>
      <td>T</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic;W/R;TGG;AGG</td>
    </tr>
  </tbody>
</table>
<p>2420 rows × 7 columns</p>
</div>




```python
likely_pathogenic_df = df_6[df_6['INFO'].str.contains("Likely_pathogenic")==True]
likely_pathogenic_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT</th>
      <th>QUAL</th>
      <th>FILTER</th>
      <th>INFO</th>
    </tr>
    <tr>
      <th>#CHROM</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>16</th>
      <td>70276973</td>
      <td>549673_1</td>
      <td>T</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;F/L;TTT;TTA</td>
    </tr>
    <tr>
      <th>16</th>
      <td>70276973</td>
      <td>549673_2</td>
      <td>T</td>
      <td>G</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;F/L;TTT;TTG</td>
    </tr>
    <tr>
      <th>9</th>
      <td>104785411</td>
      <td>1323811_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;Q/H;CAG;CAT</td>
    </tr>
    <tr>
      <th>16</th>
      <td>2286744</td>
      <td>932907_1</td>
      <td>C</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;N/K;AAC;AAA</td>
    </tr>
    <tr>
      <th>16</th>
      <td>2288150</td>
      <td>932908_1</td>
      <td>G</td>
      <td>T</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;L/F;TTG;TTT</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>12</th>
      <td>32755717</td>
      <td>1056_2</td>
      <td>T</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic/Likely_pathogenic;F/L;TTC;CTC</td>
    </tr>
    <tr>
      <th>7</th>
      <td>76329934</td>
      <td>438806_1</td>
      <td>C</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;D/E;GAC;GAA</td>
    </tr>
    <tr>
      <th>2</th>
      <td>144389856</td>
      <td>290241_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;C/S;TGC;TCC</td>
    </tr>
    <tr>
      <th>13</th>
      <td>99985328</td>
      <td>203372_1</td>
      <td>T</td>
      <td>A</td>
      <td>.</td>
      <td>.</td>
      <td>Likely_pathogenic;H/Q;CAT;CAA</td>
    </tr>
    <tr>
      <th>20</th>
      <td>45950360</td>
      <td>998091_1</td>
      <td>G</td>
      <td>C</td>
      <td>.</td>
      <td>.</td>
      <td>Pathogenic/Likely_pathogenic;G/R;GGG;CGG</td>
    </tr>
  </tbody>
</table>
<p>2750 rows × 7 columns</p>
</div>



Сохраняем наши любимые датафреймы в VCF


```python
save_VCF(benign_df, "benign.vcf")
save_VCF(likely_benign_df, "likely_benign.vcf")
save_VCF(pathogenic_df, "pathogenic.vcf")
save_VCF(likely_pathogenic_df, "likely_pathogenic.vcf")
```

## 5. Описательная статистика полученных результатов 

Посчитаем количество клинически значимых групп CLNSIG (Benign, Likely_benign, Pathogenic, Likely_pathogenic)


```python
df_1.rename(columns={'#Uploaded_variation': 'ID'}, inplace=True)
```


```python
for_search_df = df_clinvar_2.merge(df_1, on='ID', how='right')
for_search_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>CHROM</th>
      <th>POS</th>
      <th>ID</th>
      <th>REF</th>
      <th>ALT_1</th>
      <th>CLNSIG</th>
      <th>Location</th>
      <th>Allele</th>
      <th>Gene</th>
      <th>Feature</th>
      <th>Feature_type</th>
      <th>Consequence</th>
      <th>cDNA_position</th>
      <th>CDS_position</th>
      <th>Protein_position</th>
      <th>Amino_acids</th>
      <th>Codons</th>
      <th>Existing_variation</th>
      <th>Extra</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>12</td>
      <td>9079271</td>
      <td>711716</td>
      <td>C</td>
      <td>T</td>
      <td>Likely_benign</td>
      <td>12:9079271</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3162</td>
      <td>3092</td>
      <td>1031</td>
      <td>R/Q</td>
      <td>cGa/cAa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>12</td>
      <td>9079672</td>
      <td>18171</td>
      <td>T</td>
      <td>C</td>
      <td>Benign</td>
      <td>12:9079672</td>
      <td>C</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3068</td>
      <td>2998</td>
      <td>1000</td>
      <td>I/V</td>
      <td>Atc/Gtc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>12</td>
      <td>9079755</td>
      <td>18172</td>
      <td>C</td>
      <td>T</td>
      <td>Benign</td>
      <td>12:9079755</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>2985</td>
      <td>2915</td>
      <td>972</td>
      <td>C/Y</td>
      <td>tGt/tAt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>12</td>
      <td>9094987</td>
      <td>18174</td>
      <td>C</td>
      <td>T</td>
      <td>Benign</td>
      <td>12:9094987</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>2181</td>
      <td>2111</td>
      <td>704</td>
      <td>R/H</td>
      <td>cGt/cAt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>12</td>
      <td>9107574</td>
      <td>729811</td>
      <td>C</td>
      <td>T</td>
      <td>Benign</td>
      <td>12:9107574</td>
      <td>T</td>
      <td>2</td>
      <td>NM_000014.6</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>899</td>
      <td>829</td>
      <td>277</td>
      <td>D/N</td>
      <td>Gac/Aac</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=-1;CANONICAL=YES;GIVEN_...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>432375</th>
      <td>5</td>
      <td>61544156</td>
      <td>155772</td>
      <td>C</td>
      <td>T</td>
      <td>Pathogenic/Likely_pathogenic</td>
      <td>5:61544156</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3502</td>
      <td>3487</td>
      <td>1163</td>
      <td>R/W</td>
      <td>Cgg/Tgg</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432376</th>
      <td>5</td>
      <td>61544160</td>
      <td>450504</td>
      <td>A</td>
      <td>G</td>
      <td>Likely_pathogenic</td>
      <td>5:61544160</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3506</td>
      <td>3491</td>
      <td>1164</td>
      <td>H/R</td>
      <td>cAc/cGc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432377</th>
      <td>5</td>
      <td>61544165</td>
      <td>1424334</td>
      <td>A</td>
      <td>G</td>
      <td>Uncertain_significance</td>
      <td>5:61544165</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3511</td>
      <td>3496</td>
      <td>1166</td>
      <td>S/G</td>
      <td>Agt/Ggt</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432378</th>
      <td>5</td>
      <td>61544201</td>
      <td>1417067</td>
      <td>A</td>
      <td>G</td>
      <td>Uncertain_significance</td>
      <td>5:61544201</td>
      <td>G</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3547</td>
      <td>3532</td>
      <td>1178</td>
      <td>T/A</td>
      <td>Acc/Gcc</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
    <tr>
      <th>432379</th>
      <td>5</td>
      <td>61544238</td>
      <td>1353194</td>
      <td>C</td>
      <td>T</td>
      <td>Uncertain_significance</td>
      <td>5:61544238</td>
      <td>T</td>
      <td>57688</td>
      <td>NM_020928.2</td>
      <td>Transcript</td>
      <td>missense_variant</td>
      <td>3584</td>
      <td>3569</td>
      <td>1190</td>
      <td>T/I</td>
      <td>aCa/aTa</td>
      <td>-</td>
      <td>IMPACT=MODERATE;STRAND=1;CANONICAL=YES;GIVEN_R...</td>
    </tr>
  </tbody>
</table>
<p>432380 rows × 19 columns</p>
</div>




```python
benign_num = for_search_df[(for_search_df['CLNSIG'].str.contains("Benign")==True) & (for_search_df['CLNSIG'].str.contains("Likely_benign")==False)].shape[0]
likely_benign_num = for_search_df[for_search_df['CLNSIG'].str.contains("Likely_benign")==True].shape[0]
pathogenic_num = for_search_df[(for_search_df['CLNSIG'].str.contains("Pathogenic")==True) & (for_search_df['CLNSIG'].str.contains("Likely_pathogenic")==False)].shape[0]
likely_pathogenic_num = for_search_df[for_search_df['CLNSIG'].str.contains("Likely_pathogenic")==True].shape[0]
print(benign_num, likely_benign_num, pathogenic_num, likely_pathogenic_num)
```

    11069 17951 16030 19196
    


```python
benign_new_num  = benign_df.shape[0]
likely_benign_new_num  = likely_benign_df.shape[0]
pathogenic_new_num  = pathogenic_df.shape[0]
likely_pathogenic_new_num  = likely_pathogenic_df.shape[0]
```


```python
statistics_1_df = pd.DataFrame()
statistics_1_df["Index"] = ['benign', 'likely_benign', 'pathogenic', 'likely_pathogenic']
statistics_1_df["Parent"] = [benign_num, likely_benign_num, pathogenic_num, likely_pathogenic_num]
statistics_1_df["New"] = [benign_new_num, likely_benign_new_num, pathogenic_new_num, likely_pathogenic_new_num]
statistics_1_df.set_index('Index', inplace=True)
statistics_1_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Parent</th>
      <th>New</th>
    </tr>
    <tr>
      <th>Index</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>benign</th>
      <td>11069</td>
      <td>1510</td>
    </tr>
    <tr>
      <th>likely_benign</th>
      <td>17951</td>
      <td>2367</td>
    </tr>
    <tr>
      <th>pathogenic</th>
      <td>16030</td>
      <td>2420</td>
    </tr>
    <tr>
      <th>likely_pathogenic</th>
      <td>19196</td>
      <td>2750</td>
    </tr>
  </tbody>
</table>
</div>



Строим график описания выборок


```python
statistics_1_df.plot.bar()
```




    <AxesSubplot:xlabel='Index'>




    
![png](output_75_1.png)
    



```python
a = len(omim_list_2) # Список OMIM генов, любезно предоставленный Михаилом Скобловым
b = search_omim_df.shape[0] # OMIM содержащиеся в clinvar_vep_ref.vcf
c = omim_missence_canonical_df.shape[0] # OMIM missence and canonical variants
d = df_6.shape[0] # OMIM после скрипта
e = df_1.shape[0] # OMIM до скрипта
f = df_5.shape[0] # до разбивания
```

OMIM missence and canonical variants


```python
statistics_2_df = pd.DataFrame()
statistics_2_df["index"] = ['OMIM list', 'OMIM clinvar ann', 'filtred', 'after script']
statistics_2_df["number"] = [a, b, c, d]
statistics_2_df.set_index('index', inplace=True)
statistics_2_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>number</th>
    </tr>
    <tr>
      <th>index</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>OMIM list</th>
      <td>3950</td>
    </tr>
    <tr>
      <th>OMIM clinvar ann</th>
      <td>1503664</td>
    </tr>
    <tr>
      <th>filtred</th>
      <td>432380</td>
    </tr>
    <tr>
      <th>after script</th>
      <td>63383</td>
    </tr>
  </tbody>
</table>
</div>




```python
statistics_2_df.plot.bar()
```




    <AxesSubplot:xlabel='index'>




    
![png](output_79_1.png)
    


График распределение клинически значимых вариантов выборки Parent


```python
other_p = e - (benign_num + likely_benign_num + pathogenic_num + likely_pathogenic_num)
```


```python
vals = [benign_num, likely_benign_num, pathogenic_num, likely_pathogenic_num, other_p]
labels = ["benign", "likely_benign", "pathogenic", "likely_pathogenic", "other"]
fig, ax = plt.subplots(figsize=(19, 8))
ax.pie(vals, labels=labels, autopct='%1.0f%%')
ax.axis("equal")
```




    (-1.1062935875570166,
     1.1002997097113847,
     -1.1066261572262222,
     1.1034816789483628)




    
![png](output_82_1.png)
    


График распределение клинически значимых вариантов после скрипта (только новые варианты)


```python
other_new = d - (benign_new_num + likely_benign_new_num + pathogenic_new_num + likely_pathogenic_new_num)
```


```python
vals = [benign_new_num, likely_benign_new_num, pathogenic_new_num, likely_pathogenic_new_num, other_new]
labels = ["benign_new", "likely_benign_new", "pathogenic_new", "likely_pathogenic_new", "other_new"]
fig, ax = plt.subplots(figsize=(19, 8))
ax.pie(vals, labels=labels, autopct='%1.0f%%')
ax.axis("equal")
```




    (-1.1066212568485116,
     1.1003153103568541,
     -1.1066249409028193,
     1.1003940820790181)




    
![png](output_85_1.png)
    


График распределение клинически значимых вариантов выборки Parent + new


```python
vals = [benign_num + benign_new_num, likely_benign_num + likely_benign_new_num, pathogenic_num + pathogenic_new_num, likely_pathogenic_num + likely_pathogenic_new_num, other_p + other_new]
labels = ["benign_sum", "likely_benign_sum", "pathogenic_sum", "likely_pathogenic_sum", "other_sum"]
fig, ax = plt.subplots(figsize=(19, 8))
ax.pie(vals, labels=labels, autopct='%1.0f%%')
ax.axis("equal")
```




    (-1.1063625235199483,
     1.100302986906274,
     -1.106634811951224,
     1.1031306521092694)




    
![png](output_87_1.png)
    


Строим диаграмму книлической значимости полученных вариантов


```python
CLNSIG_list = df_clinvar.CLNSIG.unique()
```


```python
def make_clnsig_df(cl_list, df):
    clinsig_list = []
    for clinsig in cl_list:
        try:
            str_list = []
            row_num = df[df['INFO'].str.contains(clinsig)].shape[0]
            str_list.append(clinsig)
            str_list.append(row_num)
        except TypeError:
            pass
        clinsig_list.append(str_list)
    return pd.DataFrame(clinsig_list)
```


```python
clinsig_df = make_clnsig_df(CLNSIG_list, df_6)
clinsig_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>0</th>
      <th>1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Uncertain_significance</td>
      <td>50573.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Likely_benign</td>
      <td>2367.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Benign</td>
      <td>2205.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Conflicting_interpretations_of_pathogenicity</td>
      <td>2724.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Pathogenic</td>
      <td>2910.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>65</th>
      <td>Conflicting_interpretations_of_pathogenicity|_...</td>
      <td>2724.0</td>
    </tr>
    <tr>
      <th>66</th>
      <td>association|_risk_factor</td>
      <td>9.0</td>
    </tr>
    <tr>
      <th>67</th>
      <td>Benign|_confers_sensitivity</td>
      <td>2205.0</td>
    </tr>
    <tr>
      <th>68</th>
      <td>Benign|_association|_confers_sensitivity</td>
      <td>2205.0</td>
    </tr>
    <tr>
      <th>69</th>
      <td>Likely_pathogenic|_association</td>
      <td>2750.0</td>
    </tr>
  </tbody>
</table>
<p>70 rows × 2 columns</p>
</div>




```python
clinsig_df.groupby([0]).sum().plot(
    kind='pie', y=1, autopct='%1.0f%%', fontsize=10, figsize=(10,10))
```




    <AxesSubplot:ylabel='1'>




    
![png](output_92_1.png)
    


Процент новых альтернативных кодонов OMIM


```python
df_6.shape[0]/df_1.shape[0] * 100
```




    14.659096165410057



Среднее количество новых альтернативных кодонов на одну аминокислоту OMIM


```python
df_6.shape[0]/df_5.shape[0]
```




    1.1770069265194705



## clinsig boxplot


```python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
```


```python
benign = 'OMIM_vep_gnomad/benign_vep.vcf'
pathogenic = 'OMIM_vep_gnomad/pathogenic_vep.vcf'
likely_benign = 'OMIM_vep_gnomad/likely_benign_vep.vcf'
likely_pathogenic = 'OMIM_vep_gnomad/likely_pathogenic_vep.vcf'
```


```python
benign = 'OMIM_vep_gnomad/benign_vep.vcf'
pathogenic = 'OMIM_vep_gnomad/pathogenic_vep.vcf'
likely_benign = 'OMIM_vep_gnomad/likely_benign_vep.vcf'
likely_pathogenic = 'OMIM_vep_gnomad/likely_pathogenic_vep.vcf'
benign = 'OMIM_vep_gnomad/benign_vep.vcf'
pathogenic = 'OMIM_vep_gnomad/pathogenic_vep.vcf'
likely_benign = 'OMIM_vep_gnomad/likely_benign_vep.vcf'
likely_pathogenic = 'OMIM_vep_gnomad/likely_pathogenic_vep.vcf'
def make_header(file_name): # функция создания заголовка
    header_list = []
    with open(file_name, "r") as file:
        while True:
            line = file.readline()
            if not line:
                break
            if line.split()[0] == "#Uploaded_variation": # сохраним заголовок
                header_list.append(line.split())
                return header_list[0]
```


```python
def make_gnomAD_AF_df(name): # функция подготовки датафрейма
    df_benign = pd.read_table(name, header=None, skiprows=[i for i in range(0,54)])
    df_benign.columns = make_header(name)
    df_2 = df_benign.iloc[:, [0,1,2,3,13]]
    gnomAD_AF = df_2[df_2['Extra'].str.contains("gnomAD_AF")==True]
    return gnomAD_AF
```


```python
def make_freq_df(name): # функция создания freq_df
    df = make_gnomAD_AF_df(name)
    freq_list = []
    for row in df.itertuples(): 
        freq_list.append(float((row[5].split(';')[5].split('=')[1])))
    freq_df = pd.DataFrame(freq_list, columns=[name.split('/')[1].split('.')[0]])
    return freq_df
```


```python
benign_freq_df = make_freq_df(benign)
likely_benign_freq_df = make_freq_df(likely_benign)
pathogenic_freq_df = make_freq_df(pathogenic)
likely_pathogenic_freq_df = make_freq_df(likely_pathogenic)
```

Убираем нулевые


```python
benign_freq_non_zero_df = benign_freq_df.loc[benign_freq_df['benign_vep'] != 0]
likely_benign_non_zero_df = likely_benign_freq_df.loc[likely_benign_freq_df['likely_benign_vep'] != 0]
pathogenic_non_zero_df = pathogenic_freq_df.loc[pathogenic_freq_df['pathogenic_vep'] != 0]
likely_pathogenic_non_zero_df = likely_pathogenic_freq_df.loc[likely_pathogenic_freq_df['likely_pathogenic_vep'] != 0]
```

Соединяем датафреймы


```python
gnomAD_freq_df = pd.concat([benign_freq_non_zero_df, likely_benign_non_zero_df, pathogenic_non_zero_df, likely_pathogenic_non_zero_df], axis=1)
gnomAD_freq_df.describe()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>benign_vep</th>
      <th>likely_benign_vep</th>
      <th>pathogenic_vep</th>
      <th>likely_pathogenic_vep</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>count</th>
      <td>116.000000</td>
      <td>275.000000</td>
      <td>94.000000</td>
      <td>100.000000</td>
    </tr>
    <tr>
      <th>mean</th>
      <td>0.035397</td>
      <td>0.052555</td>
      <td>0.004591</td>
      <td>0.000109</td>
    </tr>
    <tr>
      <th>std</th>
      <td>0.133527</td>
      <td>0.127765</td>
      <td>0.007252</td>
      <td>0.000196</td>
    </tr>
    <tr>
      <th>min</th>
      <td>0.000004</td>
      <td>0.000008</td>
      <td>0.000029</td>
      <td>0.000029</td>
    </tr>
    <tr>
      <th>25%</th>
      <td>0.000055</td>
      <td>0.000030</td>
      <td>0.000042</td>
      <td>0.000029</td>
    </tr>
    <tr>
      <th>50%</th>
      <td>0.000089</td>
      <td>0.000192</td>
      <td>0.000087</td>
      <td>0.000064</td>
    </tr>
    <tr>
      <th>75%</th>
      <td>0.001440</td>
      <td>0.039880</td>
      <td>0.013380</td>
      <td>0.000069</td>
    </tr>
    <tr>
      <th>max</th>
      <td>0.618000</td>
      <td>0.708700</td>
      <td>0.023290</td>
      <td>0.001731</td>
    </tr>
  </tbody>
</table>
</div>




```python
fig, ax= plt.subplots(figsize=(10, 10))
ax.set_yscale('log')
boxplot = sns.boxplot(x="variable", y="value", data=pd.melt(gnomAD_freq_df))
boxplot.axes.set_title("gnomAD_frequency", fontsize=16)
boxplot.set_xlabel("ClinSIG_group", fontsize=14)
boxplot.set_ylabel("gnomAD_AF", fontsize=14)
plt.show()
```


    
![png](output_108_0.png)
    


Считаем процент для ненулевых clinsig


```python
b_per = benign_freq_non_zero_df.shape[0]/benign_freq_df.shape[0]*100
lb_per = likely_benign_non_zero_df.shape[0]/likely_benign_freq_df.shape[0]*100
p_per = pathogenic_non_zero_df.shape[0]/pathogenic_freq_df.shape[0]*100
lp_per = likely_pathogenic_non_zero_df.shape[0]/likely_pathogenic_freq_df.shape[0]*100
```

Считаем абсолютное количество


```python
b_abs = benign_freq_df.shape[0]
lb_abs = likely_benign_freq_df.shape[0]
p_abs = pathogenic_freq_df.shape[0]
lp_abs = likely_pathogenic_freq_df.shape[0]
```


```python
pd.DataFrame({'clinsig':['benign', 'likely_benign', 'pathogenic', 'likely_pathogenic'],
        'abs_all':[b_abs, lb_abs, p_abs, lp_abs],
        'non_zero_percent':[b_per, lb_per, p_per, lp_per]})
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>clinsig</th>
      <th>abs_all</th>
      <th>non_zero_percent</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>benign</td>
      <td>815</td>
      <td>14.233129</td>
    </tr>
    <tr>
      <th>1</th>
      <td>likely_benign</td>
      <td>1609</td>
      <td>17.091361</td>
    </tr>
    <tr>
      <th>2</th>
      <td>pathogenic</td>
      <td>876</td>
      <td>10.730594</td>
    </tr>
    <tr>
      <th>3</th>
      <td>likely_pathogenic</td>
      <td>994</td>
      <td>10.060362</td>
    </tr>
  </tbody>
</table>
</div>




```python

```
