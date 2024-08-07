# Extracting RMS values from the provided text
text = """
RMS calculated for file 1 of 28
RMS for file 1topborosil.10um.s4.tb4.8.000.xyz is: 0.224112735377 (nm)
0.0468598133958
RMS calculated for file 2 of 28
RMS for file 1topborosil.2um.s4.tb4.8.001.xyz is: 0.216471276145 (nm)
0.0618162141597
RMS calculated for file 3 of 28
RMS for file 2botborosil.10um.s2.tb4.4.024.xyz is: 0.248628667212 (nm)
0.0523662753032
RMS calculated for file 4 of 28
RMS for file 2botborosil.2um.s1.tb4.4.023.xyz is: 0.228836787478 (nm)
0.0601167810115
RMS calculated for file 5 of 28
RMS for file 2botborosil.2um.s2.tb4.4.025.xyz is: 0.245187236641 (nm)
0.0644464943852
RMS calculated for file 6 of 28
RMS for file 2midborosil.2um.s1.tb4.4.019.xyz is: 0.253863141053 (nm)
0.0594905827514
RMS calculated for file 7 of 28
RMS for file 2midborosil.2um.s2.tb4.4.021.xyz is: 0.243906914111 (nm)
0.056225479606
RMS calculated for file 8 of 28
RMS for file 2midborosilica.10um.s4.tb4.8.008.difference.xyz is: 0.237119125348 (nm)
0.056225479606
RMS calculated for file 9 of 28
RMS for file 2midborosilica.10um.s4.tb4.8.008.xyz is: 0.237119125348 (nm)
0.0563747274849
RMS calculated for file 10 of 28
RMS for file 2midborosilica.10um.s4.tb4.8.009.xyz is: 0.237433627536 (nm)
0.0526588870443
RMS calculated for file 11 of 28
RMS for file 2topborosil.10um.s4.tb4.8.002.xyz is: 0.22947524277 (nm)
0.0431574103649
RMS calculated for file 12 of 28
RMS for file 2topborosil.2um.s1.tb4.2.014.xyz is: 0.207743616905 (nm)
0.0518204225079
RMS calculated for file 13 of 28
RMS for file 2topborosil.2um.s2.tb4.2.017.xyz is: 0.227640994788 (nm)
0.0616026322434
RMS calculated for file 14 of 28
RMS for file 2topborosil.2um.s3.tb4.4.028.xyz is: 0.248198775669 (nm)
0.0475154877159
RMS calculated for file 15 of 28
RMS for file 2topborosil.2um.s4.tb4.8.003.xyz is: 0.217980475538 (nm)
0.0586677858688
RMS calculated for file 16 of 28
RMS for file borosil.10um.s1.tb4.1.000.xyz is: 0.242214338694 (nm)
0.0650847133937
RMS calculated for file 17 of 28
RMS for file botborosil.10um.s1.tb4.1.007.xyz is: 0.255117058218 (nm)
0.0653392800667
RMS calculated for file 18 of 28
RMS for file botborosil.10um.s2.tb4.2.011.xyz is: 0.255615492619 (nm)
0.0626837181563
RMS calculated for file 19 of 28
RMS for file botborosil.2um.s1.tb4.1.008.xyz is: 0.25036716669 (nm)
0.060010338198
RMS calculated for file 20 of 28
RMS for file botborosil.2um.s2.tb4.2.012.xyz is: 0.244970076128 (nm)
0.07351031789
RMS calculated for file 21 of 28
RMS for file midborosil.10um.s1.tb4.1.004.xyz is: 0.271127862622 (nm)
0.0718805463758
RMS calculated for file 22 of 28
RMS for file midborosil.10um.s2.tb4.1.006.xyz is: 0.268105476214 (nm)
0.0540986482672
RMS calculated for file 23 of 28
RMS for file midborosil.10um.s3.tb4.2.009.xyz is: 0.232591161197 (nm)
0.0629380120732
RMS calculated for file 24 of 28
RMS for file midborosil.2um.s1.tb4.1.005.xyz is: 0.250874494665 (nm)
0.0472393283959
RMS calculated for file 25 of 28
RMS for file midborosil.2um.s3.tb4.2.010.xyz is: 0.217346102785 (nm)
0.0550230797413
RMS calculated for file 26 of 28
RMS for file topborosil.10um.s2.tb4.1.002.xyz is: 0.234569989004 (nm)
0.0601504279256
RMS calculated for file 27 of 28
RMS for file topborosil.2um.s1.tb4.1.001.xyz is: 0.245255841777 (nm)
0.0561820159983
RMS calculated for file 28 of 28
RMS for file topborosil.2um.s2.tb4.1.003.xyz is: 0.237027458321 (nm)
Average RMS: 0.239603580745
"""

# Extracting all the RMS values
import re

# Regular expressions to extract RMS values for 2um and 10um files separately
rms_2um = re.findall(r"RMS for file .+2um.+ is: ([0-9.]+) \(nm\)", text)
rms_10um = re.findall(r"RMS for file .+10um.+ is: ([0-9.]+) \(nm\)", text)

# Converting the extracted values to floats
rms_2um_float = [float(value) for value in rms_2um]
rms_10um_float = [float(value) for value in rms_10um]

# Print the separated RMS values
print("RMS values for 2um files:", rms_2um_float)
print("RMS values for 10um files:", rms_10um_float)

import numpy as np

# Calculate the average RMS values
print("Average RMS for 2um files:", np.mean(rms_2um_float))
print("Average RMS for 10um files:", np.mean(rms_10um_float))

# Calculate the standard deviation of the RMS values
print("Standard deviation of RMS for 2um files:", np.std(rms_2um_float))
print("Standard deviation of RMS for 10um files:", np.std(rms_10um_float))
