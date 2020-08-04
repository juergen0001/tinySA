/*
 * RDA5815.c - RDA5815 library for Arduino
 *
 * Copyright (C) 2015 - 2016 Jason Milldrum <milldrum@gmail.com>
 *                           Dana H. Myers <k6jq@comcast.net>
 *
 * Some tuning algorithms derived from clk-RDA5815.c in the Linux kernel.
 * Sebastian Hesselbarth <sebastian.hesselbarth@gmail.com>
 * Rabeeh Khoury <rabeeh@solid-run.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hal.h"
#include "nanovna.h"

#define RDA5815_BUS_BASE_ADDR            12
 
static int RDA5815_write(uint8_t reg, uint8_t data);

void  RDA5815_init(void)
{
    my_microsecond_delay(1000);
// Chip register soft reset
RDA5815_write(0x04,0x04);
RDA5815_write(0x04,0x05);
        my_microsecond_delay(1000);
//pll setting
RDA5815_write(0x1a,0x13);
RDA5815_write(0x38,0x90);
RDA5815_write(0x39,0x15);
RDA5815_write(0x3A,0x00);
RDA5815_write(0x3B,0x00);
RDA5815_write(0x0c,0xE2);
RDA5815_write(0x2e,0x6F);
RDA5815_write(0x72,0x07);
RDA5815_write(0x73,0x60);
RDA5815_write(0x74,0x76);

RDA5815_write(0x5b,0x20);
RDA5815_write(0x2f,0x57);
RDA5815_write(0x0d,0x70);
RDA5815_write(0x16,0x03);
RDA5815_write(0x18,0x4B);
RDA5815_write(0x30,0xFF);
RDA5815_write(0x5c,0xFF);
RDA5815_write(0x6c,0xFF);
RDA5815_write(0x6e,0xFF);
RDA5815_write(0x65,0x80);
RDA5815_write(0x70,0x3F);
RDA5815_write(0x71,0x3F);
RDA5815_write(0x75,0x06);
RDA5815_write(0x76,0x40);
RDA5815_write(0x77,0x89);
RDA5815_write(0x53,0xA8);
RDA5815_write(0x46,0x21);
RDA5815_write(0x47,0x84);
RDA5815_write(0x48,0x10);
RDA5815_write(0x49,0x08);
RDA5815_write(0x60,0x80);
RDA5815_write(0x61,0x80);
RDA5815_write(0x6A,0x08);
RDA5815_write(0x6B,0x63);
RDA5815_write(0x69,0xF8);
RDA5815_write(0x57,0x74);
RDA5815_write(0x05,0x88);
RDA5815_write(0x06,0x88);
RDA5815_write(0x15,0xAA);
RDA5815_write(0x4a,0xbb);
RDA5815_write(0x4b,0xbb);
//agc setting
RDA5815_write(0x4f,0x40);
RDA5815_write(0x5b,0x20);
// for blocker
RDA5815_write(0x16,0x10);//stage setting
RDA5815_write(0x18,0x20);
RDA5815_write(0x30,0x30);
RDA5815_write(0x5c,0x30);
RDA5815_write(0x6c,0x30);
RDA5815_write(0x6e,0x70);
RDA5815_write(0x1b,0xB2);
RDA5815_write(0x1d,0xB2);
RDA5815_write(0x1f,0xB2);
RDA5815_write(0x21,0xB2);
RDA5815_write(0x23,0xB6);
RDA5815_write(0x25,0xB6);
RDA5815_write(0x27,0xBA);
RDA5815_write(0x29,0xBF);
RDA5815_write(0xb3,0xFF);
RDA5815_write(0xb5,0xFF);
RDA5815_write(0x17,0xF0);
RDA5815_write(0x19,0xF0);
RDA5815_write(0x31,0xF0);
RDA5815_write(0x5d,0xF1);
RDA5815_write(0x6d,0xF2);
RDA5815_write(0x6f,0xF2);
RDA5815_write(0x1c,0x31);
RDA5815_write(0x1e,0x72);
RDA5815_write(0x20,0x96);
RDA5815_write(0x22,0xBA);
RDA5815_write(0x24,0xBA);
RDA5815_write(0x26,0xBE);
RDA5815_write(0x28,0xCE);
RDA5815_write(0x2a,0xDE);
RDA5815_write(0xb4,0x0F);
RDA5815_write(0xb6,0x0F);
RDA5815_write(0xb7,0x10); //start
RDA5815_write(0xb9,0x10);
RDA5815_write(0xbb,0x00);
RDA5815_write(0xbd,0x00);
RDA5815_write(0xbf,0x00);
RDA5815_write(0xc1,0x10);
RDA5815_write(0xc3,0x10);
RDA5815_write(0xc5,0x10);
RDA5815_write(0xa3,0x19);
RDA5815_write(0xa5,0x2E);
RDA5815_write(0xa7,0x37);
RDA5815_write(0xa9,0x47);
RDA5815_write(0xab,0x47);
RDA5815_write(0xad,0x3F);
RDA5815_write(0xaf,0x00);
RDA5815_write(0xb1,0x37);
RDA5815_write(0xb8,0x47); //end
RDA5815_write(0xba,0x3F);
RDA5815_write(0xbc,0x37);
RDA5815_write(0xbe,0x3F);
RDA5815_write(0xc0,0x3F);
RDA5815_write(0xc2,0x3F);
RDA5815_write(0xc4,0x3F);
RDA5815_write(0xc6,0x3F);
RDA5815_write(0xa4,0x47);
RDA5815_write(0xa6,0x57);
RDA5815_write(0xa8,0x5F);
RDA5815_write(0xaa,0x70);
RDA5815_write(0xac,0x70);
RDA5815_write(0xae,0x68);
RDA5815_write(0xb0,0x37);
RDA5815_write(0xb2,0x37);
RDA5815_write(0x81,0x77); //rise
RDA5815_write(0x82,0x68);
RDA5815_write(0x83,0x70);
RDA5815_write(0x84,0x68);
RDA5815_write(0x85,0x68);
RDA5815_write(0x86,0x68);
RDA5815_write(0x87,0x70);
RDA5815_write(0x88,0x47);
RDA5815_write(0x89,0x68);
RDA5815_write(0x8a,0x8E);
RDA5815_write(0x8b,0x8E);
RDA5815_write(0x8c,0x8E);
RDA5815_write(0x8d,0x9C);
RDA5815_write(0x8e,0x9A);
RDA5815_write(0x8f,0x37);
RDA5815_write(0x90,0x00); //fall
RDA5815_write(0x91,0x00);
RDA5815_write(0x92,0x00);
RDA5815_write(0x93,0x00);
RDA5815_write(0x94,0x00);
RDA5815_write(0x95,0x00);
RDA5815_write(0x96,0x00);
RDA5815_write(0x97,0x00);
RDA5815_write(0x98,0x00);
RDA5815_write(0x99,0x00);
RDA5815_write(0x9a,0x10);
RDA5815_write(0x9b,0x24);
RDA5815_write(0x9c,0x10);
RDA5815_write(0x9d,0x00);
RDA5815_write(0x9e,0x00);

#if 1        
// 30MHz crystal
RDA5815_write(0x72,0x06);    // 0x0 6 A
RDA5815_write(0x73,0xA0);
RDA5815_write(0x74,0x6A);

RDA5815_write(0x75,0x05);      // 0x0 5 A
RDA5815_write(0x76,0xA0);     // 0x0 7 B
RDA5815_write(0x77,0x7B);

RDA5815_write(0x79,0x04);   // 0x04 7A AA AA
RDA5815_write(0x7A,0x7A);
RDA5815_write(0x7B,0xAA);
RDA5815_write(0x7C,0xAA);
#endif        
        
        my_microsecond_delay(10000); // Initial configuration ready}

RDA5815_write(0x04,0xc1);
RDA5815_write(0x07,0x07);
RDA5815_write(0x08,0x68);
RDA5815_write(0x09,0x4b);
RDA5815_write(0x0a,0xda);
RDA5815_write(0x0b,0x0a);
RDA5815_write(0x04,0xc3);

  my_microsecond_delay(5000);//Wait 5ms;


}


uint8_t RDA5815_set_freq(uint32_t fPLL, uint32_t fSym )
{			                   
 	uint8_t buffer;
 	uint32_t temp_value = 0;
	uint32_t bw;/*,temp_value1 = 0,temp_value2=0 ;*/
	uint8_t Filter_bw_control_bit;	

	static uint32_t old_fPLL = 0;
	if (old_fPLL == fPLL)
	  return 1;
	old_fPLL = fPLL;

	fPLL -= 30;             // Subtract 30kHz for XCO error
	
    uint32_t MHz = (fPLL / 1000);
	uint32_t kHz = fPLL - MHz*1000;

	RDA5815_write(0x04,0xc1); //add by rda 2011.8.9,RXON = 0 , change normal working state to idle state
  //RDA5815_write(0x2b,0x95);//clk_interface_27m=0  add by rda 2012.1.12     
	//set frequency start
//	temp_value = (uint32_t)fPLL * /* 77672; // */ 69905;//((2<<21) / RDA5815_XTALFREQ);
    temp_value = MHz * 69905 + kHz * 70;

	buffer = ((uint8_t)((temp_value>>24)&0xff));
	RDA5815_write(0x07,buffer);
	buffer = ((uint8_t)((temp_value>>16)&0xff));	
	RDA5815_write(0x08,buffer);	
   	buffer = ((uint8_t)((temp_value>>8)&0xff));
	RDA5815_write(0x09,buffer);	
   	buffer = ((uint8_t)( temp_value&0xff));
	RDA5815_write(0x0a,buffer);
	//set frequency end
	
	// set Filter bandwidth start
	bw=fSym;
	
	if(bw<4000)
	bw= 4000;    // KHz
	else if(bw>45000)
	bw = 40000;   // KHz    //modify by rda 2012.1.12   
	
	Filter_bw_control_bit = (uint8_t)((bw*135/200+4000)/1000);
	
	Filter_bw_control_bit&=0x3f;
	RDA5815_write(0x0b,Filter_bw_control_bit);
	// set Filter bandwidth end
	
	RDA5815_write(0x04,0xc3); //add by rda 2011.8.9,RXON = 0 ,rxon=1,normal working
	//RDA5815_write(0x2b,0x97);//clk_interface_27m=1  add by rda 2012.1.12  
	my_microsecond_delay(5000);//Wait 5ms;

	return 1;   
}

static int RDA5815_write(uint8_t reg, uint8_t data)
{
  int addr = RDA5815_BUS_BASE_ADDR;
  uint8_t buf[] = { reg, data };
  i2cAcquireBus(&I2CD1);
  msg_t mr = i2cMasterTransmitTimeout(&I2CD1, addr, buf, 2, NULL, 0, 1000);
  i2cReleaseBus(&I2CD1);
  return mr == MSG_OK;
}
