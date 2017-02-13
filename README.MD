#Introduction

+	Spatt includes some 'Large Deviation' samples.  Migrated this over to win64
to play with those samples
+	Included latest gsl at this time.

#Library folder's structure

* `\gsl-2.1-src` - 
* `\gsl-2.1-src\VC14` -   (open libgsl.sln)
* `\spatt-1.2.2.t3` -      (open spatt-1.2.2.t3.sln)

#Build Instructions

+	Build libgsl.sln first. Review that configuration is for a x64 build. This creates a new directory with required include files
\gsl-2.1-src\gsl ( this subdirectory is created during this build)
+	for gsl-2.1-src, I did not run full regression tests
+	Only test was whether it worked with spatt-1.2.2.t3
+	Build spatt-1.2.2.t3.sln
+	Used VC14 (Visual Studio 2015)
+	VC14, Update 1 works fine with relative paths.Previous release was buggy.

#Changes to original code. Spatt

+	In spatt-1.2.2.t3\src  are several file names, renamed with prefix 'main'.
	Each of these have a int main function.  These are test files.
	In VC14, Solution Explorer, you can include or exclude as to test.
	
#Changes to original code. gsl

See file "Changes.to.gsl-2.1.to.build.in.VC14.64.bat.txt"

#Output from:
hiv.fasta -l 4 -m 1 -V --all-words --max-pvalue 0.1 > test_result_ldspatt.txt

count length  			4			
Markov order 1 specified						
length("hiv.fasta")			10901		includes white spaces	
count memory usage (bytes) 			1360			
sequence length			9718		w/o white spaces	
						
						
pattern	nobs	natt	stat			
aaac	42	59.66	-1.314358			
aaat	96	75.83	1.14361			
aacc	26	42.62	-1.677783			
aagg	64	93.14	-2.368286			
aatt	83	58.88	1.998653			
agaa	148	115.55	1.851555			
aggc	33	55.08	-2.319139			
aggt	25	52.06	-3.910779			
agta	96	59	4.271916			
agtc	17	29.5	-1.382496			
agtt	33	47.27	-1.081144			
atcc	38	26.25	1.022044			
atga	50	67.27	-1.066825			
atgc	22	37.49	-1.669311			
caaa	65	84.52	-1.129944			
ccaa	40	60.38	-1.766736			
cccg	14	5.15	2.246298			
cgta	1	5.47	-1.207353			
ctaa	30	48.44	-1.821535			
ctgt	45	23.35	3.495926			
ctta	24	37.61	-1.25949			
gaca	78	60.48	1.055798			
gagg	47	68.66	-1.655312			
gagt	31	44.88	-1.075787			
gcac	22	33.71	-1.034677			
ggaa	103	74.56	2.226494			
ggcg	14	6.83	1.249208			
gtac	42	22.46	3.006728			
gtag	59	41	1.499082			
gtcc	5	15.9	-2.245011			
gtga	23	40.74	-2.041267			
gtgg	50	32.84	1.665666			
gttc	9	17.29	-1.060902			
gttg	18	29.82	-1.18263			
taag	36	66.93	-3.903585			
tagg	42	57.25	-1.012561			
tagt	64	37.41	3.340595			
tcat	19	32.04	-1.343134			
tctg	40	24.7	1.769076			
tctt	12	22.94	-1.368687			
tgaa	33	66.52	-4.699309			
tgag	35	61.26	-3.0141			
tgga	84	56.89	2.564963			
tggg	62	45.86	1.154112			
tgtg	48	29.3	1.987367			
ttga	28	51.66	-2.919791			
