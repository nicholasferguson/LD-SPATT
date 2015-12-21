(1)	Migrated LD-SPATT and dependencies over to win64
	Included latest gsl at this time.
	
(2) Used VC14 (Visual Studio 2015)
	It is a bit buggy with relative paths in a project's properties.
(3) In spatt-1.2.2.t3\src  are several file names prefixed by 'main'.
	Each of these have a int main function.
	In VC14, Solution Explorer, you can include or exclude to test.
(4) only two solution files were tested
		spatt-1.2.2.t3\spatt-1.2.2.t3.sln		
		gsl-2.1-src\VC14\libgsl.sln
(5) for gsl-2.1-src, I did not run full regression tests.

	
LD-SPATT does include a "large deviation"