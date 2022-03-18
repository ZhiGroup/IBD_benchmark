IBD Benchmark project source code.

List of file:

Program.cs         Program entrance, load configurations, check inputs, then run program.
GroupCaller.cs     Group calling module, organized to call core functions to run a full set of test.
BinCalculator.cs   Core calculation module, detailed calculation methords could be found here.
Loader.cs          Loading unit, parsers could be modified here for differet formats.
Utl.cs             Supporting utility module.



Compiling and Execution

Precompiled versions for ubuntu 18, ubuntu 21 and windows 10 are already uploaded to directory ```/IBD_Benchmark_ExperiencePackage/```. However if those version do not work for you, you can compile the source code. As following:

For windows:

Download code file, compile with Microsoft Visual Studio, use Release version for better performance.

For Linux and macOS:

Dowload the code files, use mono compiler ```mcs -out:IBD_BM *.cs```
If you have not installed mono compiler, please follow the guide in https://www.mono-project.com/ to install.

Once compiled, you may run:
```./IBD_BM``` for linux
```mono IBD_BM``` for macOS

Note, input file paths and most configurations could be modified in Config.txt, it may not be necessary to modify in the code.
