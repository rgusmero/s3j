# s3j
### c routines for 3j symbol evaluation

The main subroutine

_**double** s3j(**double** j1, **double** j2, **double** j3, **double** m1, **double** m2, **double** m3)_

calculates the following symbol:
```
    ( j1 j2 j3 )
   (            ) = delta(m1+m2+m3,0) * (-1)^(j1-j2-m3) * 
    ( m1 m2 m3 )

      +-
      |  (j1+j2-j3)! (j1-j2+j3)! (-j1+j2+j3)! 
    * | -------------------------------------- ...
      |
      +-
                                                                 -+ 1/2
           (j1-m1)! (j1+m1)! (j2-m2)! (j2+m2)! (j3-m3)! (j3+m3)!  |
      ... ------------------------------------------------------- |     * 
                              (j1+j2+j3+1)!                       |
                                                                 -+

         +---
          \                       (-1)^k
      *    |   ---------------------------------------------------------------------
          /      k! (j1+j2-j3-k)! (j1-m1-k)! (j2+m2-k)! (j3-j2+m1+k)! (j3-j1-m2+k)!
         +---
           k
```           
### Configuration
The maximum allowed factorial is S3J_MAX_FACT (currently 25!).

### Testing
Defining the macro S3J_TEST enables the compilation of a very small test suite.

### License

**Author: Riccardo Gusmeroli, Ph.D. (rgusmero@elet.polimi.it)**

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. 
