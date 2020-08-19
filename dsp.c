/*
 * Copyright (c) 2014-2015, TAKAHASHI Tomohiro (TTRFTECH) edy555@gmail.com
 * All rights reserved.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * The software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include <arm_math.h>
#include "nanovna.h"

#ifdef __DUMP_CMD__
int16_t samp_buf[SAMPLE_LEN];
int16_t ref_buf[SAMPLE_LEN];
#endif //__DUMP_CMD__

#ifdef __FLOAT_FFT__
void FFT(float data[], int m, bool forward);
float data[MAX_FFT_LEN*2];
float window[MAX_FFT_LEN/2];
#endif

#ifdef __INT_FFT__
static int fix_fft(short fr[], short fi[], short m, short inverse);
int16_t data[MAX_FFT_LEN*2];
int16_t window[MAX_FFT_LEN/2];
#endif

int fft_len_div_2;
int fft_len;
int fft_fill = 0;
int fft_bits = 1;
static int16_t *l_c;

#define SHOW(f,i)   {char string_buf[12];plot_printf(string_buf, sizeof string_buf, "%f", f);ili9341_drawstring(string_buf, 150 + i*80  , FREQUENCIES_YPOS);}

#define CAVER 2
#define complex_mul(r,a,b) { r[0] = (a)[0]*(b)[0] - (a)[1]*(b)[1]; r[1] = (a)[0]*(b)[1] + (a)[1]*(b)[1]; }
// (a + bi) / (c + di) = ((ac + bd) / (c2 + d2)) + ((bc - ad) / (c2 + d2)i
#define complex_div(r,a,b) { float cd2 = (b)[0]*(b)[0] + (b)[1]*(b)[1]; if (cd2>0) { (r)[0] = ((a)[0]*(b)[0] + (a)[1]*(b)[1])/cd2; (r)[1] = ((a)[1]*(b)[0] - (a)[0]*(b)[1])/cd2; } }

static float phase(float *v)
{
  return 2 * atan2f(v[1], v[0]) / 3.141592653 * 90;
}

static float ampl(float *t)
{
  return sqrt(t[0]*t[0] + t[1]*t[1]);
}

float last_max = 0;
float current_max = 0;

volatile int32_t ph1,ph2,ph3;
#define LENGTH ((1<<15) - 1)
volatile float c1=0,c2=0.97;
//int32_t aizero, arzero;
float aizero, arzero;

void calculate_correlation(void)
{
  current_max = 0;
  ph1 = 0;
  ph2 = 0;
  ph3 = 0;

  arzero = 0;
  aizero = 0;
  for (int i = 0; i < fft_len; i++)
  {
    int k = i*2;
    aizero += data[k];
    arzero += data[k+1];
  }
  aizero /= fft_len;
  arzero /= fft_len;

  for (int i = 0; i < fft_len; i++)
  {
    int k = i*2;
    float s_i = data[k] - aizero,
          s_q = data[k+1] - arzero;

//    ph1 +=  - sign(i) * s_q;
//    ph2 += sign(i) * s_i;
    if (s_i > 0) {
      ph1 -= s_q;
      ph2 += s_i;
    } else {
      ph1 += s_q;
      ph2 -= s_i;
    }
//    ph3 += sign(q) * s_q;
    if (s_q > 0)
      ph3 += s_q;
    else
      ph3 -= s_q;

    int16_t w = i;       // calculate index in window
    if (w >= fft_len_div_2)
      w = fft_len - 1 - w;

    data[k] = s_i * window[w] * c2;     // Window and I/Q balance compensate
    data[k+1] = s_q * window[w];
  }

  int aver;
  if (ph2 > 100000) {
    if (c2 == 1.0) aver = 0; else aver = 40;
    c1 = ((c1 * aver) + (float)ph1/(float)ph2 ) / (aver + 1);
    if (c2 < 0.9) c2 = 0.9;
    if (c2 > 1.1) c2 = 1.1;
    c2 = ((c2 * aver) + sqrt(((float)ph3*ph3 - (float)ph1*ph1)/((float)ph2*ph2)) ) / (aver+1);
//    SHOW(c1,0);
//    SHOW(c2,1);
  }
#ifdef __FLOAT_FFT__                                            // real FFT
  FFT(data, fft_bits, true);
#endif

#ifdef __INT_FFT__
  int16_t *re = &data[0];
  int16_t *im = &data[length/2];

  for (int i = 0, j = 0; i < (int)length; i +=2, j++) {
    re[j] = (((int32_t)(capture[i] - rzero)) * window[j]) >> 16;  //
    im[j] = (((int32_t)(capture[i+1] - izero)) * window[j]) >> 16;  //
  }
  int scale = fix_fft(re,im, bits,false);
#endif


}






volatile float rsum = 0, isum = 0;
//volatile int32_t s_max;

void dsp_process(int16_t *capture, size_t length)   // 2 samples per fft point
{
  int len_div_2 = length/2;
//  int32_t rzero = 0, izero = 0;

  l_c = capture;            // Only to get raw audio

//  s_max = 0;
  for (int i=0;i<len_div_2;i++) {
    int16_t s;
    s = capture[i*2+0];
//    rzero += s;
//    s -= arzero;
    data[fft_fill++] = (float)s;
//    if (s < 0)
//      s = -s;
//    if (s_max < s)
//      s_max = s;
    s = capture[i*2+1];
//    izero += s;
//    s -= aizero;
    data[fft_fill++] = (float)s;
//    if (s < 0)
//      s = -s;
//    if (s_max < s)
//      s_max = s;
  }
  if (fft_fill < fft_len*2)
    wait_count++;           // Get one more buffer
  else
    fft_fill = 0;

#define AVER_ZERO  100
//  arzero = (arzero * AVER_ZERO * len_div_2 + rzero)/(AVER_ZERO+1)/len_div_2;
//  aizero = (aizero * AVER_ZERO * len_div_2 + izero)/(AVER_ZERO+1)/len_div_2;
}



static int dsp_filled = false;
static int dsp_index = 0;

int audio_level = 0;

void set_audio_level(int l)
{
  audio_level = l;
}

void dsp_fill(void)
{
  wait_count = 3;
  while (wait_count) __WFI();

  if (audio_level==0)
    calculate_correlation();
  dsp_filled = true;
  dsp_index = 0;
}

float dsp_get(int fi)
{
  int samples = ((fft_len) / 256);
  if (samples < 1)
    samples = 1;

  if (fft_len > 256)           // subsample
    fi = fi * samples;

  if (fi > fft_len)            // no more data
    return -150.0;

  if (fi < fft_len/2)
    fi = fi + fft_len/2;
  else
    fi = fi - fft_len/2 + 1;
//  fi = fft_len - fi;     // Invert frequencies due to I/Q phase error

  float sub_data,max_data = 0;
  while (samples --) {
#ifdef __FLOAT_FFT__
    sub_data = data[2*fi]*data[2*fi] + data[2*fi+1]*data[2*fi+1];         // dBm
#endif
#ifdef __INT_FFT__
  sub_data = data[fi]*data[fi] + data[fi+fft_len]*data[fi+fft_len];         // dBm
#endif
    if (max_data < sub_data)
      max_data = sub_data;
    fi++;
  }
  return max_data;
}

float dsp_get_one(int fi)
{
  if (fi > fft_len || fi < 0)            // no more data
    return -150.0;

  if (fi < fft_len/2)
    fi = fi + fft_len/2;
  else
    fi = fi - fft_len/2 + 1;
//  fi = fft_len - fi;     // Invert frequencies due to I/Q phase error

#ifdef __FLOAT_FFT__
  return data[2*fi]*data[2*fi] + data[2*fi+1]*data[2*fi+1];         // dBm
#endif
#ifdef __INT_FFT__
  return data[fi]*data[fi] + data[fi+fft_len]*data[fi+fft_len];         // dBm
#endif
}


float dsp_getmax(int fft_steps, int fft_step) {
  float maxdata = 0;
  if (dsp_filled) {
    if (dsp_index == fft_len - 1)
      dsp_filled = false;
    if (audio_level)
      return data[2*(dsp_index++) + (audio_level & 1)] / (audio_level/2); // Get the raw audio signal left or right channel
    maxdata = dsp_get(dsp_index++);
  } else {
    if (fft_step == 0) {
      wait_count = 3;
      while (wait_count) __WFI();
      calculate_correlation();
    }
    float submax;
    int fft_bucket = (fft_len + fft_steps - 1) / fft_steps;
    int i = -fft_bucket/2;
    do {
      submax = dsp_get_one((fft_step * fft_len) / fft_steps + i);
      if (maxdata < submax)
        maxdata = submax;
      i++;
    } while (i < fft_bucket/2 );
  }
#if 1
#else
  calculate_correlation();
#define IGNORE  10
  for (int i=IGNORE; i < fft_len - IGNORE; i++) {
    float sub_data = dsp_get(i);
    if (maxdata < sub_data)
      maxdata = sub_data;
    }
  }
#endif
#ifdef __FLOAT_FFT__
    float RSSI = 10*log10(maxdata/(fft_len*fft_len)) - 10;         // dBm
#endif
#ifdef __INT_FFT__
    float RSSI = 10*log10(maxdata) - 90;         // dBm
#endif
  return RSSI;
}

void dsp_init(int len) {
  if (len < MIN_FFT_LEN)
    len = MIN_FFT_LEN;
  if (len > MAX_FFT_LEN)
    len = MAX_FFT_LEN;
  if (fft_len == len)
    return;
  fft_len = len;
  for (int i = 0; i < fft_len/2; i++) {
//#define PI  3.14159265358979
#define A 0.16
#define WINDOW(n) ((1.0-A)/2 - 0.5 * cos((2.0*PI*n)/(fft_len-1)) + (A/2.0)*cos((4*PI*n)/(fft_len-1)))
#ifdef __FLOAT_FFT__
    window[i] = WINDOW(i);
#else
    window[i] = (int16_t) (MAX_INT16 * WINDOW(i));
#endif
  }
  fft_len_div_2 = fft_len/2;
  fft_bits = 1;
  while (1<<fft_bits < fft_len)
    fft_bits++;

}

#ifdef __FLOAT_FFT__           // float fft

void FFT(float data[], int m, bool forward)
{
    int n, i, i1, j, k, i2, l, l1, l2;
    float c1, c2, tx, ty, t1, t2, u1, u2, z;

    n = 1<<m;
//  n /= 2;

    // Do the bit reversal
    i2 = n >> 1;
    j = 0;
    for (i = 0; i < n - 1; i++)
    {
        if (i < j)
        {
            tx = data[i * 2];
            ty = data[i * 2 + 1];
            data[i * 2] = data[j * 2];
            data[i * 2 + 1] = data[j * 2 + 1];
            data[j * 2] = tx;
            data[j * 2 + 1] = ty;
        }
        k = i2;

        while (k <= j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    // Compute the FFT
    c1 = -1.0f;
    c2 = 0.0f;
    l2 = 1;
    for (l = 0; l < m; l++)
    {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0f;
        u2 = 0.0f;
        for (j = 0; j < l1; j++)
        {
            for (i = j; i < n; i += l2)
            {
                i1 = i + l1;
                t1 = u1 * data[i1 * 2] - u2 * data[i1 * 2 + 1];
                t2 = u1 * data[i1 * 2 + 1] + u2 * data[i1 * 2];
                data[i1 * 2] = data[i * 2] - t1;
                data[i1 * 2 + 1] = data[i * 2 + 1] - t2;
                data[i * 2] += t1;
                data[i * 2 + 1] += t2;
            }
            z = u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        c2 = sqrt((1.0f - c1) / 2.0f);
        if (!forward)
            c2 = -c2;
        c1 = sqrt((1.0f + c1) / 2.0f);
    }

    // Scaling for inverse transform
    if (!forward)
    {
        for (i = 0; i < n; i++)
        {
            data[i * 2] /= n;
            data[i * 2 + 1] /= n;
        }
    }
}

#endif

#ifdef __INT_FFT__

/* fix_fft.c - Fixed-point in-place Fast Fourier Transform  */
/*
  All data are fixed-point short integers, in which -32768
  to +32768 represent -1.0 to +1.0 respectively. Integer
  arithmetic is used for speed, instead of the more natural
  floating-point.

  For the forward FFT (time -> freq), fixed scaling is
  performed to prevent arithmetic overflow, and to map a 0dB
  sine/cosine wave (i.e. amplitude = 32767) to two -6dB freq
  coefficients. The return value is always 0.

  For the inverse FFT (freq -> time), fixed scaling cannot be
  done, as two 0dB coefficients would sum to a peak amplitude
  of 64K, overflowing the 32k range of the fixed-point integers.
  Thus, the fix_fft() routine performs variable scaling, and
  returns a value which is the number of bits LEFT by which
  the output must be shifted to get the actual amplitude
  (i.e. if fix_fft() returns 3, each value of fr[] and fi[]
  must be multiplied by 8 (2**3) for proper scaling.
  Clearly, this cannot be done within fixed-point short
  integers. In practice, if the result is to be used as a
  filter, the scale_shift can usually be ignored, as the
  result will be approximately correctly normalized as is.

  Written by:  Tom Roberts  11/8/89
  Made portable:  Malcolm Slaney 12/15/94 malcolm@interval.com
  Enhanced:  Dimitrios P. Bouras  14 Jun 2006 dbouras@ieee.org
*/

#define N_WAVE      1024    /* full length of Sinewave[] */
#define LOG2_N_WAVE 10      /* log2(N_WAVE) */

/*
  Henceforth "short" implies 16-bit word. If this is not
  the case in your architecture, please replace "short"
  with a type definition which *is* a 16-bit word.
*/

/*
  Since we only use 3/4 of N_WAVE, we define only
  this many samples, in order to conserve data space.
*/
short Sinewave[N_WAVE-N_WAVE/4] = {
      0,    201,    402,    603,    804,   1005,   1206,   1406,
   1607,   1808,   2009,   2209,   2410,   2610,   2811,   3011,
   3211,   3411,   3611,   3811,   4011,   4210,   4409,   4608,
   4807,   5006,   5205,   5403,   5601,   5799,   5997,   6195,
   6392,   6589,   6786,   6982,   7179,   7375,   7571,   7766,
   7961,   8156,   8351,   8545,   8739,   8932,   9126,   9319,
   9511,   9703,   9895,  10087,  10278,  10469,  10659,  10849,
  11038,  11227,  11416,  11604,  11792,  11980,  12166,  12353,
  12539,  12724,  12909,  13094,  13278,  13462,  13645,  13827,
  14009,  14191,  14372,  14552,  14732,  14911,  15090,  15268,
  15446,  15623,  15799,  15975,  16150,  16325,  16499,  16672,
  16845,  17017,  17189,  17360,  17530,  17699,  17868,  18036,
  18204,  18371,  18537,  18702,  18867,  19031,  19194,  19357,
  19519,  19680,  19840,  20000,  20159,  20317,  20474,  20631,
  20787,  20942,  21096,  21249,  21402,  21554,  21705,  21855,
  22004,  22153,  22301,  22448,  22594,  22739,  22883,  23027,
  23169,  23311,  23452,  23592,  23731,  23869,  24006,  24143,
  24278,  24413,  24546,  24679,  24811,  24942,  25072,  25201,
  25329,  25456,  25582,  25707,  25831,  25954,  26077,  26198,
  26318,  26437,  26556,  26673,  26789,  26905,  27019,  27132,
  27244,  27355,  27466,  27575,  27683,  27790,  27896,  28001,
  28105,  28208,  28309,  28410,  28510,  28608,  28706,  28802,
  28897,  28992,  29085,  29177,  29268,  29358,  29446,  29534,
  29621,  29706,  29790,  29873,  29955,  30036,  30116,  30195,
  30272,  30349,  30424,  30498,  30571,  30643,  30713,  30783,
  30851,  30918,  30984,  31049,  31113,  31175,  31236,  31297,
  31356,  31413,  31470,  31525,  31580,  31633,  31684,  31735,
  31785,  31833,  31880,  31926,  31970,  32014,  32056,  32097,
  32137,  32176,  32213,  32249,  32284,  32318,  32350,  32382,
  32412,  32441,  32468,  32495,  32520,  32544,  32567,  32588,
  32609,  32628,  32646,  32662,  32678,  32692,  32705,  32717,
  32727,  32736,  32744,  32751,  32757,  32761,  32764,  32766,
  32767,  32766,  32764,  32761,  32757,  32751,  32744,  32736,
  32727,  32717,  32705,  32692,  32678,  32662,  32646,  32628,
  32609,  32588,  32567,  32544,  32520,  32495,  32468,  32441,
  32412,  32382,  32350,  32318,  32284,  32249,  32213,  32176,
  32137,  32097,  32056,  32014,  31970,  31926,  31880,  31833,
  31785,  31735,  31684,  31633,  31580,  31525,  31470,  31413,
  31356,  31297,  31236,  31175,  31113,  31049,  30984,  30918,
  30851,  30783,  30713,  30643,  30571,  30498,  30424,  30349,
  30272,  30195,  30116,  30036,  29955,  29873,  29790,  29706,
  29621,  29534,  29446,  29358,  29268,  29177,  29085,  28992,
  28897,  28802,  28706,  28608,  28510,  28410,  28309,  28208,
  28105,  28001,  27896,  27790,  27683,  27575,  27466,  27355,
  27244,  27132,  27019,  26905,  26789,  26673,  26556,  26437,
  26318,  26198,  26077,  25954,  25831,  25707,  25582,  25456,
  25329,  25201,  25072,  24942,  24811,  24679,  24546,  24413,
  24278,  24143,  24006,  23869,  23731,  23592,  23452,  23311,
  23169,  23027,  22883,  22739,  22594,  22448,  22301,  22153,
  22004,  21855,  21705,  21554,  21402,  21249,  21096,  20942,
  20787,  20631,  20474,  20317,  20159,  20000,  19840,  19680,
  19519,  19357,  19194,  19031,  18867,  18702,  18537,  18371,
  18204,  18036,  17868,  17699,  17530,  17360,  17189,  17017,
  16845,  16672,  16499,  16325,  16150,  15975,  15799,  15623,
  15446,  15268,  15090,  14911,  14732,  14552,  14372,  14191,
  14009,  13827,  13645,  13462,  13278,  13094,  12909,  12724,
  12539,  12353,  12166,  11980,  11792,  11604,  11416,  11227,
  11038,  10849,  10659,  10469,  10278,  10087,   9895,   9703,
   9511,   9319,   9126,   8932,   8739,   8545,   8351,   8156,
   7961,   7766,   7571,   7375,   7179,   6982,   6786,   6589,
   6392,   6195,   5997,   5799,   5601,   5403,   5205,   5006,
   4807,   4608,   4409,   4210,   4011,   3811,   3611,   3411,
   3211,   3011,   2811,   2610,   2410,   2209,   2009,   1808,
   1607,   1406,   1206,   1005,    804,    603,    402,    201,
      0,   -201,   -402,   -603,   -804,  -1005,  -1206,  -1406,
  -1607,  -1808,  -2009,  -2209,  -2410,  -2610,  -2811,  -3011,
  -3211,  -3411,  -3611,  -3811,  -4011,  -4210,  -4409,  -4608,
  -4807,  -5006,  -5205,  -5403,  -5601,  -5799,  -5997,  -6195,
  -6392,  -6589,  -6786,  -6982,  -7179,  -7375,  -7571,  -7766,
  -7961,  -8156,  -8351,  -8545,  -8739,  -8932,  -9126,  -9319,
  -9511,  -9703,  -9895, -10087, -10278, -10469, -10659, -10849,
 -11038, -11227, -11416, -11604, -11792, -11980, -12166, -12353,
 -12539, -12724, -12909, -13094, -13278, -13462, -13645, -13827,
 -14009, -14191, -14372, -14552, -14732, -14911, -15090, -15268,
 -15446, -15623, -15799, -15975, -16150, -16325, -16499, -16672,
 -16845, -17017, -17189, -17360, -17530, -17699, -17868, -18036,
 -18204, -18371, -18537, -18702, -18867, -19031, -19194, -19357,
 -19519, -19680, -19840, -20000, -20159, -20317, -20474, -20631,
 -20787, -20942, -21096, -21249, -21402, -21554, -21705, -21855,
 -22004, -22153, -22301, -22448, -22594, -22739, -22883, -23027,
 -23169, -23311, -23452, -23592, -23731, -23869, -24006, -24143,
 -24278, -24413, -24546, -24679, -24811, -24942, -25072, -25201,
 -25329, -25456, -25582, -25707, -25831, -25954, -26077, -26198,
 -26318, -26437, -26556, -26673, -26789, -26905, -27019, -27132,
 -27244, -27355, -27466, -27575, -27683, -27790, -27896, -28001,
 -28105, -28208, -28309, -28410, -28510, -28608, -28706, -28802,
 -28897, -28992, -29085, -29177, -29268, -29358, -29446, -29534,
 -29621, -29706, -29790, -29873, -29955, -30036, -30116, -30195,
 -30272, -30349, -30424, -30498, -30571, -30643, -30713, -30783,
 -30851, -30918, -30984, -31049, -31113, -31175, -31236, -31297,
 -31356, -31413, -31470, -31525, -31580, -31633, -31684, -31735,
 -31785, -31833, -31880, -31926, -31970, -32014, -32056, -32097,
 -32137, -32176, -32213, -32249, -32284, -32318, -32350, -32382,
 -32412, -32441, -32468, -32495, -32520, -32544, -32567, -32588,
 -32609, -32628, -32646, -32662, -32678, -32692, -32705, -32717,
 -32727, -32736, -32744, -32751, -32757, -32761, -32764, -32766,
};

/*
  FIX_MPY() - fixed-point multiplication & scaling.
  Substitute inline assembly for hardware-specific
  optimization suited to a particluar DSP processor.
  Scaling ensures that result remains 16-bit.
*/
inline short FIX_MPY(short a, short b)
{
    /* shift right one less bit (i.e. 15-1) */
    int c = ((int)a * (int)b) >> 14;
    /* last bit shifted out = rounding-bit */
    b = c & 0x01;
    /* last shift + rounding bit */
    a = (c >> 1) + b;
    return a;
}

/*
  fix_fft() - perform forward/inverse fast Fourier transform.
  fr[n],fi[n] are real and imaginary arrays, both INPUT AND
  RESULT (in-place FFT), with 0 <= n < 2**m; set inverse to
  0 for forward transform (FFT), or 1 for iFFT.
*/
int fix_fft(short fr[], short fi[], short m, short inverse)
{
    int mr, nn, i, j, l, k, istep, n, scale, shift;
    short qr, qi, tr, ti, wr, wi;

    n = 1 << m;

    /* max FFT size = N_WAVE */
    if (n > N_WAVE)
        return -1;

    mr = 0;
    nn = n - 1;
    scale = 0;

    /* decimation in time - re-order data */
    for (m=1; m<=nn; ++m) {
        l = n;
        do {
            l >>= 1;
        } while (mr+l > nn);
        mr = (mr & (l-1)) + l;

        if (mr <= m)
            continue;
        tr = fr[m];
        fr[m] = fr[mr];
        fr[mr] = tr;
        ti = fi[m];
        fi[m] = fi[mr];
        fi[mr] = ti;
    }

    l = 1;
    k = LOG2_N_WAVE-1;
    while (l < n) {
        if (inverse) {
            /* variable scaling, depending upon data */
            shift = 0;
            for (i=0; i<n; ++i) {
                j = fr[i];
                if (j < 0)
                    j = -j;
                m = fi[i];
                if (m < 0)
                    m = -m;
                if (j > 16383 || m > 16383) {
                    shift = 1;
                    break;
                }
            }
            if (shift)
                ++scale;
        } else {
            /*
              fixed scaling, for proper normalization --
              there will be log2(n) passes, so this results
              in an overall factor of 1/n, distributed to
              maximize arithmetic accuracy.
            */
            shift = 1;
        }
        /*
          it may not be obvious, but the shift will be
          performed on each data point exactly once,
          during this pass.
        */
        istep = l << 1;
        for (m=0; m<l; ++m) {
            j = m << k;
            /* 0 <= j < N_WAVE/2 */
            wr =  Sinewave[j+N_WAVE/4];
            wi = -Sinewave[j];
            if (inverse)
                wi = -wi;
            if (shift) {
                wr >>= 1;
                wi >>= 1;
            }
            for (i=m; i<n; i+=istep) {
                j = i + l;
                tr = FIX_MPY(wr,fr[j]) - FIX_MPY(wi,fi[j]);
                ti = FIX_MPY(wr,fi[j]) + FIX_MPY(wi,fr[j]);
                qr = fr[i];
                qi = fi[i];
                if (shift) {
                    qr >>= 1;
                    qi >>= 1;
                }
                fr[j] = qr - tr;
                fi[j] = qi - ti;
                fr[i] = qr + tr;
                fi[i] = qi + ti;
            }
        }
        --k;
        l = istep;
    }
    return scale;
}

/*
  fix_fftr() - forward/inverse FFT on array of real numbers.
  Real FFT/iFFT using half-size complex FFT by distributing
  even/odd samples into real/imaginary arrays respectively.
  In order to save data space (i.e. to avoid two arrays, one
  for real, one for imaginary samples), we proceed in the
  following two steps: a) samples are rearranged in the real
  array so that all even samples are in places 0-(N/2-1) and
  all imaginary samples in places (N/2)-(N-1), and b) fix_fft
  is called with fr and fi pointing to index 0 and index N/2
  respectively in the original array. The above guarantees
  that fix_fft "sees" consecutive real samples as alternating
  real and imaginary samples in the complex array.
*/
int fix_fftr(short f[], int m, int inverse)
{
    int i, N = 1<<(m-1), scale = 0;
    short tt, *fr=f, *fi=&f[N];

    if (inverse)
        scale = fix_fft(fi, fr, m-1, inverse);
    for (i=1; i<N; i+=2) {
        tt = f[N+i-1];
        f[N+i-1] = f[i];
        f[i] = tt;
    }
    if (! inverse)
        scale = fix_fft(fi, fr, m-1, inverse);
    return scale;
}

#endif
