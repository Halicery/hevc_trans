/*************************************************************************************
*   Fast 1-D Inverse Integer Transform for HEVC
* 
* 
*   Reduced number of multiplications by recursive 
*   factorization and circulant transform of the odd-part. 
*   Generated code based on the paper
*   http://halicery.com/paper/Fast%20algorithm%20for%20the%201-D%20integer%20transform%20of%20HEVC%20with%20circulants.pdf
*   using the following Permutation Vectors: 
* 
* 		2,5,14,-9,7,-13,6,-16,15,-12,3,8,-10,4,11,-1,
* 		8,7,4,6,1,-2,5,3,
* 		3,4,2,-1,
* 
* 
*   Copyright (c) 2017 Attila Tarpai 
*   
*   Permission is hereby granted, free of charge, to any person obtaining a copy
*   of this software and associated documentation files (the "Software"), to deal
*   in the Software without restriction, including without limitation the rights
*   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*   copies of the Software, and to permit persons to whom the Software is
*   furnished to do so, subject to the following conditions:
*   
*   The above copyright notice and this permission notice shall be included in all
*   copies or substantial portions of the Software.
*   
*   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*   SOFTWARE.
*/



////////  ODD-MATRIX CIRCULANT SOLVERS  //////////////

static void r2(int *x, int *y)
{
 y[1]= (x[0] + x[1]) * 36;
 y[0]= y[1] + x[0] * 47;
 y[1]+= x[1] * -119;
}

static void r4(int *x, int *y)
{
y[0]= x[0]+x[2];
y[1]= x[1]+x[3];
 y[3]= (y[0] + y[1]) * -50;
 y[2]= y[3] + y[0] * -39;
 y[3]+= y[1] * 32;
 y[1]= (x[0] + x[1]) * 125;
 y[0]= y[2] + y[1] + x[0] * -18;
 y[1]= y[3] + y[1] + x[1] * -196;
 x[1]= (x[2] + x[3]) * -25;
 y[2]+= x[1] + x[2] * 96;
 y[3]+= x[1] + x[3] * 132;
}

static void r8(int *x, int *y)
{
y[0]= x[0]+x[4];
y[1]= x[1]+x[5];
y[2]= x[2]+x[6];
y[3]= x[3]+x[7];
y[4]= y[0]+y[2];
y[5]= y[1]+y[3];
 y[7]= (y[4] + y[5]) * 43;
 y[6]= y[7] + y[4] * 27;
 y[7]+= y[5] * 47;
 y[5]= (y[0] + y[1]) * -18;
 y[4]= y[6] + y[5] + y[0] * -43;
 y[5]= y[7] + y[5] + y[1] * -2;
 y[1]= (y[2] + y[3]) * -130;
 y[6]+= y[1] + y[2] * 150;
 y[7]+= y[1] + y[3] * 97;
y[0]= x[0]+x[2];
y[1]= x[1]+x[3];
 y[3]= (y[0] + y[1]) * -123;
 y[2]= y[3] + y[0] * -4;
 y[3]+= y[1] * 42;
 y[1]= (x[0] + x[1]) * 185;
 y[0]= y[4] + y[2] + y[1] + x[0] * -157;
 y[1]= y[5] + y[3] + y[1] + x[1] * -231;
 x[1]= (x[2] + x[3]) * 235;
 y[2]= y[6] + y[2] + x[1] + x[2] * -189;
 y[3]= y[7] + y[3] + x[1] + x[3] * -141;
x[0]= x[4]+x[6];
x[1]= x[5]+x[7];
 x[3]= (x[0] + x[1]) * 37;
 x[2]= x[3] + x[0] * -50;
 x[3]+= x[1] * -136;
 x[1]= (x[4] + x[5]) * -149;
 y[4]+= x[2] + x[1] + x[4] * 243;
 y[5]+= x[3] + x[1] + x[5] * 235;
 x[5]= (x[6] + x[7]) * 25;
 y[6]+= x[2] + x[5] + x[6] * -111;
 y[7]+= x[3] + x[5] + x[7] * -53;
}

static void r16(int *x, int *y)
{
y[0]= x[0]+x[8];
y[1]= x[1]+x[9];
y[2]= x[2]+x[10];
y[3]= x[3]+x[11];
y[4]= x[4]+x[12];
y[5]= x[5]+x[13];
y[6]= x[6]+x[14];
y[7]= x[7]+x[15];
y[8]= y[0]+y[4];
y[9]= y[1]+y[5];
y[10]= y[2]+y[6];
y[11]= y[3]+y[7];
y[12]= y[8]+y[10];
y[13]= y[9]+y[11];
 y[15]= (y[12] + y[13]) * -90;
 y[14]= y[15] + y[12] * 0;
 y[15]+= y[13] * 8;
 y[13]= (y[8] + y[9]) * 136;
 y[12]= y[14] + y[13] + y[8] * 39;
 y[13]= y[15] + y[13] + y[9] * -144;
 y[9]= (y[10] + y[11]) * 68;
 y[14]+= y[9] + y[10] * -60;
 y[15]+= y[9] + y[11] * 75;
y[8]= y[0]+y[2];
y[9]= y[1]+y[3];
 y[11]= (y[8] + y[9]) * 36;
 y[10]= y[11] + y[8] * 121;
 y[11]+= y[9] * 131;
 y[9]= (y[0] + y[1]) * 6;
 y[8]= y[12] + y[10] + y[9] + y[0] * -286;
 y[9]= y[13] + y[11] + y[9] + y[1] * -16;
 y[1]= (y[2] + y[3]) * 32;
 y[10]= y[14] + y[10] + y[1] + y[2] * -22;
 y[11]= y[15] + y[11] + y[1] + y[3] * -350;
y[0]= y[4]+y[6];
y[1]= y[5]+y[7];
 y[3]= (y[0] + y[1]) * 17;
 y[2]= y[3] + y[0] * 134;
 y[3]+= y[1] * 96;
 y[1]= (y[4] + y[5]) * -85;
 y[12]+= y[2] + y[1] + y[4] * -233;
 y[13]+= y[3] + y[1] + y[5] * 123;
 y[5]= (y[6] + y[7]) * -73;
 y[14]+= y[2] + y[5] + y[6] * 35;
 y[15]+= y[3] + y[5] + y[7] * -97;
y[0]= x[0]+x[4];
y[1]= x[1]+x[5];
y[2]= x[2]+x[6];
y[3]= x[3]+x[7];
y[4]= y[0]+y[2];
y[5]= y[1]+y[3];
 y[7]= (y[4] + y[5]) * 103;
 y[6]= y[7] + y[4] * -17;
 y[7]+= y[5] * -59;
 y[5]= (y[0] + y[1]) * -71;
 y[4]= y[6] + y[5] + y[0] * -131;
 y[5]= y[7] + y[5] + y[1] * 113;
 y[1]= (y[2] + y[3]) * 7;
 y[6]+= y[1] + y[2] * -49;
 y[7]+= y[1] + y[3] * -45;
y[0]= x[0]+x[2];
y[1]= x[1]+x[3];
 y[3]= (y[0] + y[1]) * 24;
 y[2]= y[3] + y[0] * -238;
 y[3]+= y[1] * -184;
 y[1]= (x[0] + x[1]) * -122;
 y[0]= y[8] + y[4] + y[2] + y[1] + x[0] * 572;
 y[1]= y[9] + y[5] + y[3] + y[1] + x[1] * 68;
 x[1]= (x[2] + x[3]) * -102;
 y[2]= y[10] + y[6] + y[2] + x[1] + x[2] * 156;
 y[3]= y[11] + y[7] + y[3] + x[1] + x[3] * 342;
x[0]= x[4]+x[6];
x[1]= x[5]+x[7];
 x[3]= (x[0] + x[1]) * -84;
 x[2]= x[3] + x[0] * 4;
 x[3]+= x[1] * 94;
 x[1]= (x[4] + x[5]) * 162;
 y[4]= y[12] + y[4] + x[2] + x[1] + x[4] * 78;
 y[5]= y[13] + y[5] + x[3] + x[1] + x[5] * -252;
 x[5]= (x[6] + x[7]) * 98;
 y[6]= y[14] + y[6] + x[2] + x[5] + x[6] * -8;
 y[7]= y[15] + y[7] + x[3] + x[5] + x[7] * -208;
x[0]= x[8]+x[12];
x[1]= x[9]+x[13];
x[2]= x[10]+x[14];
x[3]= x[11]+x[15];
x[4]= x[0]+x[2];
x[5]= x[1]+x[3];
 x[7]= (x[4] + x[5]) * 77;
 x[6]= x[7] + x[4] * 17;
 x[7]+= x[5] * 43;
 x[5]= (x[0] + x[1]) * -201;
 x[4]= x[6] + x[5] + x[0] * 53;
 x[5]= x[7] + x[5] + x[1] * 175;
 x[1]= (x[2] + x[3]) * -143;
 x[6]+= x[1] + x[2] * 169;
 x[7]+= x[1] + x[3] * -105;
x[0]= x[8]+x[10];
x[1]= x[9]+x[11];
 x[3]= (x[0] + x[1]) * -96;
 x[2]= x[3] + x[0] * -4;
 x[3]+= x[1] * -78;
 x[1]= (x[8] + x[9]) * 110;
 y[8]+= x[4] + x[2] + x[1] + x[8] * 0;
 y[9]+= x[5] + x[3] + x[1] + x[9] * -36;
 x[9]= (x[10] + x[11]) * 38;
 y[10]+= x[6] + x[2] + x[9] + x[10] * -112;
 y[11]+= x[7] + x[3] + x[9] + x[11] * 358;
x[8]= x[12]+x[14];
x[9]= x[13]+x[15];
 x[11]= (x[8] + x[9]) * 50;
 x[10]= x[11] + x[8] * -272;
 x[11]+= x[9] * -286;
 x[9]= (x[12] + x[13]) * 8;
 y[12]+= x[4] + x[10] + x[9] + x[12] * 388;
 y[13]+= x[5] + x[11] + x[9] + x[13] * 6;
 x[13]= (x[14] + x[15]) * 48;
 y[14]+= x[6] + x[10] + x[13] + x[14] * -62;
 y[15]+= x[7] + x[11] + x[13] + x[15] * 402;
}



///////  SUB-MATRIX SOLVERS  //////////////

static void t4(int *x, int *y)
{
	int r[4];
	r[0]= ( x[0] + x[1] ) * 64;
	r[1]= ( x[0] - x[1] ) * 64;
	r2(x+2, r+2);
	y[0] = r[0] + r[2];
	y[1] = r[1] + r[3];
	y[2] = r[1] - r[3];
	y[3] = r[0] - r[2];
}

static void t8(int *x, int *y)
{
	int r[8];
	t4(x, r);
	r4(x+4, r+4);
	y[0] = r[0] - r[7];
	y[1] = r[1] + r[6];
	y[2] = r[2] + r[4];
	y[3] = r[3] + r[5];
	y[4] = r[3] - r[5];
	y[5] = r[2] - r[4];
	y[6] = r[1] - r[6];
	y[7] = r[0] + r[7];
}

static void t16(int *x, int *y)
{
	int r[16];
	t8(x, r);
	r8(x+8, r+8);
	y[0] = r[0] + r[12];
	y[1] = r[1] - r[13];
	y[2] = r[2] + r[15];
	y[3] = r[3] + r[10];
	y[4] = r[4] + r[14];
	y[5] = r[5] + r[11];
	y[6] = r[6] + r[9];
	y[7] = r[7] + r[8];
	y[8] = r[7] - r[8];
	y[9] = r[6] - r[9];
	y[10] = r[5] - r[11];
	y[11] = r[4] - r[14];
	y[12] = r[3] - r[10];
	y[13] = r[2] - r[15];
	y[14] = r[1] + r[13];
	y[15] = r[0] - r[12];
}

static void t32(int *x, int *y)
{
	int r[32];
	t16(x, r);
	r16(x+16, r+16);
	y[0] = r[0] - r[31];
	y[1] = r[1] + r[16];
	y[2] = r[2] + r[26];
	y[3] = r[3] + r[29];
	y[4] = r[4] + r[17];
	y[5] = r[5] + r[22];
	y[6] = r[6] + r[20];
	y[7] = r[7] + r[27];
	y[8] = r[8] - r[19];
	y[9] = r[9] - r[28];
	y[10] = r[10] + r[30];
	y[11] = r[11] - r[25];
	y[12] = r[12] - r[21];
	y[13] = r[13] + r[18];
	y[14] = r[14] + r[24];
	y[15] = r[15] - r[23];
	y[16] = r[15] + r[23];
	y[17] = r[14] - r[24];
	y[18] = r[13] - r[18];
	y[19] = r[12] + r[21];
	y[20] = r[11] + r[25];
	y[21] = r[10] - r[30];
	y[22] = r[9] + r[28];
	y[23] = r[8] + r[19];
	y[24] = r[7] - r[27];
	y[25] = r[6] - r[20];
	y[26] = r[5] - r[22];
	y[27] = r[4] - r[17];
	y[28] = r[3] - r[29];
	y[29] = r[2] - r[26];
	y[30] = r[1] - r[16];
	y[31] = r[0] + r[31];
}




//////////////  ENTRY POINTS   ///////////////////

extern void tr4(int *x, int *y)
{
	int a[4];
	a[0] = x[0];
	a[1] = x[2];
	a[2] = x[1];
	a[3] = x[3];
	t4(a, y);
}

extern void tr8(int *x, int *y)	
{
	int a[8];
	a[0] = x[0];
	a[1] = x[4];
	a[2] = x[2];
	a[3] = x[6];
	a[4] = x[5];
	a[5] = x[7];
	a[6] = x[3];
	a[7] = -x[1];
	t8(a, y);
}

extern void tr16(int *x, int *y)
{
	int a[16];
	a[0] = x[0];
	a[1] = x[8];
	a[2] = x[4];
	a[3] = x[12];
	a[4] = x[10];
	a[5] = x[14];
	a[6] = x[6];
	a[7] = -x[2];
	a[8] = x[15];
	a[9] = x[13];
	a[10] = x[7];
	a[11] = x[11];
	a[12] = x[1];
	a[13] = -x[3];
	a[14] = x[9];
	a[15] = x[5];
	t16(a, y);
}

extern void tr32(int *x, int *y)
{
	int a[32];
	a[0] = x[0];
	a[1] = x[16];
	a[2] = x[8];
	a[3] = x[24];
	a[4] = x[20];
	a[5] = x[28];
	a[6] = x[12];
	a[7] = -x[4];
	a[8] = x[30];
	a[9] = x[26];
	a[10] = x[14];
	a[11] = x[22];
	a[12] = x[2];
	a[13] = -x[6];
	a[14] = x[18];
	a[15] = x[10];
	a[16] = x[3];
	a[17] = x[9];
	a[18] = x[27];
	a[19] = -x[17];
	a[20] = x[13];
	a[21] = -x[25];
	a[22] = x[11];
	a[23] = -x[31];
	a[24] = x[29];
	a[25] = -x[23];
	a[26] = x[5];
	a[27] = x[15];
	a[28] = -x[19];
	a[29] = x[7];
	a[30] = x[21];
	a[31] = -x[1];
	t32(a, y);
}
