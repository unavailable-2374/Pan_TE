KAnalyze is a program to convert DNA and RNA sequence data into overlapping segments of the same length, k-mers. For example, 8-mers have 8 bases each.

KAnalyze requires Java 7 (also known as 1.7). To check for the correct version, run:
java -version

If you do not see version 1.7.0 or later, please update Java before continuing.


There are two modes, count and stream. Count outputs a sorted tab-delimited file of k-mers and their counts. Stream writes a list of k-mers as they are read.

The count mode has a graphical interface. To launch it, double-click on "countgui" or run it from the command line with no options. The count command also launches the graphical interface from the command line if it is started without any command line options.

To count k-mers in the included test file, test.fa, run:
count -k 8 -f fasta -o test.kc test.fa

The above command counts all 8-mers (-k 8) in the input file (test.fa) of format FASTA (-f fasta) and writes a tab-delimited file of the k-mers and their counts to test.kc (-o test.kc).

For command line help, run:
count -h

For full documentation on KAnalyze, see the included KAnalyzeManual.pdf (requires a PDF reader).


Copyright (c) 2014 Peter A. Audano III

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Library General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program; see the file COPYING.LESSER.  If not, see
<http://www.gnu.org/licenses/>
