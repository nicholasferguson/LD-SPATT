(1)	Migrated LD-SPATT and dependencies over to win64
	Included latest gsl at this time.
(1.A) To build.  First build gsl then spatt
	
(2) Used VC14 (Visual Studio 2015)
	VC14, Update 1 works fine with relative paths.Previous release was buggy.
(3) In spatt-1.2.2.t3\src  are several file names, renamed with prefix 'main'.
	Each of these have a int main function.  These are test files.
	In VC14, Solution Explorer, you can include or exclude as to test.
(4) only two solution files were used in builds.
		spatt-1.2.2.t3\spatt-1.2.2.t3.sln		
		gsl-2.1-src\VC14\libgsl.sln
(5) for gsl-2.1-src, I did not run full regression tests.
	Only test was whether it worked with spatt-1.2.2.t3

	
LD-SPATT does include a "large deviation test. 