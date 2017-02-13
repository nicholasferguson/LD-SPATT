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

#Changes to original code

+	In spatt-1.2.2.t3\src  are several file names, renamed with prefix 'main'.
	Each of these have a int main function.  These are test files.
	In VC14, Solution Explorer, you can include or exclude as to test.
	
