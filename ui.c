/* Copyright (c) 2014-2015, TAKAHASHI Tomohiro (TTRFTECH) edy555@gmail.com
 * Copyright (c) 2020-2022, Erik Kaashoek
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

#include "ch.h"
#include "hal.h"
#include "chprintf.h"
#include "nanovna.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#pragma GCC push_options
#pragma GCC optimize ("Os")

uistat_t uistat = {
 current_trace: 0,
 lever_mode: LM_MARKER,
 marker_delta: FALSE,
 marker_noise: FALSE,
 marker_tracking : FALSE,
 auto_center_marker : FALSE,
 text : "",
};

#define NO_EVENT                    0
#define EVT_BUTTON_SINGLE_CLICK     0x01
#define EVT_BUTTON_DOUBLE_CLICK     0x02
#define EVT_BUTTON_DOWN_LONG        0x04
#define EVT_UP                  0x10
#define EVT_DOWN                0x20
#define EVT_REPEAT              0x40

#define BUTTON_DOWN_LONG_TICKS      MS2ST(500)   // 500ms
#define BUTTON_DOUBLE_TICKS         MS2ST(250)   // 250ms
#define BUTTON_REPEAT_TICKS         MS2ST( 40)   //  40ms
#define BUTTON_DEBOUNCE_TICKS       MS2ST(  2)   //   2ms

/* lever switch assignment */
#define BIT_UP1     3
#define BIT_PUSH    2
#define BIT_DOWN1   1

#define READ_PORT() palReadPort(GPIOA)
#define BUTTON_MASK 0b1110

static uint16_t last_button = 0b0000;
static uint32_t last_button_down_ticks;
static uint32_t last_button_repeat_ticks;

#define MENU_USE_AUTOHEIGHT
#ifdef MENU_USE_AUTOHEIGHT
static uint16_t menu_button_height = MENU_BUTTON_HEIGHT_N(MENU_BUTTON_MIN);
#endif

volatile uint8_t operation_requested = OP_NONE;

int8_t previous_marker = MARKER_INVALID;

enum {
  UI_NORMAL, UI_MENU, UI_KEYPAD
};

#define NUMINPUT_LEN 12
static uint8_t ui_mode = UI_NORMAL;
static uint8_t keypad_mode;
static char    kp_buf[NUMINPUT_LEN+1];
static int8_t  kp_index = 0;
static char   *kp_help_text = NULL;
static uint8_t menu_current_level = 0;
static int  selection = 0;

static const uint8_t slider_bitmap[]=
{
  _BMP8(0b11111110),
  _BMP8(0b11111110),
  _BMP8(0b01111100),
  _BMP8(0b00111000),
  _BMP8(0b00010000)
};

// Button definition (used in MT_ADV_CALLBACK for custom)
#define BUTTON_ICON_NONE            -1
#define BUTTON_ICON_NOCHECK          0
#define BUTTON_ICON_CHECK            1
#define BUTTON_ICON_CHECK_AUTO       2
#define BUTTON_ICON_CHECK_MANUAL     3
#define BUTTON_ICON_GROUP            4
#define BUTTON_ICON_GROUP_CHECKED    5

#define CHECK_ICON(S) ((S) ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK)
#define GROUP_ICON(S) ((S) ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP)
#define AUTO_ICON(S)  (S>=2?BUTTON_ICON_CHECK_AUTO:S)            // Depends on order of ICONs!!!!!

#define BUTTON_BORDER_NONE           0x00
#define BUTTON_BORDER_WIDTH_MASK     0x0F

// Define mask for draw border (if 1 use light color, if 0 dark)
#define BUTTON_BORDER_TYPE_MASK      0xF0
#define BUTTON_BORDER_TOP            0x10
#define BUTTON_BORDER_BOTTOM         0x20
#define BUTTON_BORDER_LEFT           0x40
#define BUTTON_BORDER_RIGHT          0x80

#define BUTTON_BORDER_FLAT           0x00
#define BUTTON_BORDER_RISE           (BUTTON_BORDER_TOP|BUTTON_BORDER_RIGHT)
#define BUTTON_BORDER_FALLING        (BUTTON_BORDER_BOTTOM|BUTTON_BORDER_LEFT)

// Touch screen
#define EVT_TOUCH_NONE     0
#define EVT_TOUCH_DOWN     1
#define EVT_TOUCH_PRESSED  2
#define EVT_TOUCH_RELEASED 3

#define TOUCH_INTERRUPT_ENABLED   1
static uint8_t touch_status_flag = 0;
static int8_t last_touch_status = EVT_TOUCH_NONE;
static int16_t last_touch_x;
static int16_t last_touch_y;

#define KP_CONTINUE 0
#define KP_DONE 1
#define KP_CANCEL 2

static void ui_mode_keypad(int _keypad_mode);
// static void draw_menu(void);
static void leave_ui_mode(void);
static void erase_menu_buttons(void);
static void ui_process_keypad(void);
static void choose_active_marker(void);
static void menu_move_back(bool leave_ui);
//static const menuitem_t menu_marker_type[];

static int btn_check(void)
{
  systime_t ticks;
  // Debounce input
  while(TRUE){
    ticks = chVTGetSystemTimeX();
    if(ticks - last_button_down_ticks > BUTTON_DEBOUNCE_TICKS)
      break;
    chThdSleepMilliseconds(10);
  }
  int status = 0;
  uint16_t cur_button = READ_PORT() & BUTTON_MASK;
  // Detect only changed and pressed buttons
  uint16_t button_set = (last_button ^ cur_button) & cur_button;
  last_button_down_ticks = ticks;
  last_button = cur_button;

  if (button_set & (1<<BIT_PUSH))
    status |= EVT_BUTTON_SINGLE_CLICK;
  if (button_set & (1<<BIT_UP1))
    status |= EVT_UP;
  if (button_set & (1<<BIT_DOWN1))
    status |= EVT_DOWN;
  return status;
}

static int btn_wait_release(void)
{
  while (TRUE) {
    systime_t ticks = chVTGetSystemTimeX();
    systime_t dt = ticks - last_button_down_ticks;
    // Debounce input
//    if (dt < BUTTON_DEBOUNCE_TICKS){
//      chThdSleepMilliseconds(10);
//      continue;
//    }
    chThdSleepMilliseconds(1);
    uint16_t cur_button = READ_PORT() & BUTTON_MASK;
    uint16_t changed = last_button ^ cur_button;
    if (dt >= BUTTON_DOWN_LONG_TICKS && (cur_button & (1<<BIT_PUSH)))
      return EVT_BUTTON_DOWN_LONG;
    else if (changed & (1<<BIT_PUSH)) { // release
      last_button = cur_button;
      last_button_down_ticks = ticks;
      return EVT_BUTTON_SINGLE_CLICK;
    }
    if (changed) {
      // finished
      last_button = cur_button;
      last_button_down_ticks = ticks;
      return 0;
    }

    if (dt > BUTTON_DOWN_LONG_TICKS &&
      ticks > last_button_repeat_ticks) {
      int status = 0;
      if (cur_button & (1<<BIT_DOWN1))
        status |= EVT_DOWN | EVT_REPEAT;
      if (cur_button & (1<<BIT_UP1))
        status |= EVT_UP | EVT_REPEAT;
      last_button_repeat_ticks = ticks + BUTTON_REPEAT_TICKS;
      return status;
    }
  }
}


#define SOFTWARE_TOUCH
//*******************************************************************************
// Software Touch module
//*******************************************************************************
#ifdef SOFTWARE_TOUCH
// ADC read count for measure X and Y (2^N count)
#define TOUCH_X_N 3
#define TOUCH_Y_N 3
static int
touch_measure_y(void)
{
  // drive low to high on X line (At this state after touch_prepare_sense)
//  palSetPadMode(GPIOB, GPIOB_XN, PAL_MODE_OUTPUT_PUSHPULL); //
//  palSetPadMode(GPIOA, GPIOA_XP, PAL_MODE_OUTPUT_PUSHPULL); //
  // drive low to high on X line (coordinates from top to bottom)
  palClearPad(GPIOB, GPIOB_XN);
//  palSetPad(GPIOA, GPIOA_XP);

  // open Y line (At this state after touch_prepare_sense)
//  palSetPadMode(GPIOB, GPIOB_YN, PAL_MODE_INPUT);        // Hi-z mode
  palSetPadMode(GPIOA, GPIOA_YP, PAL_MODE_INPUT_ANALOG);   // <- ADC_TOUCH_Y channel

//  chThdSleepMilliseconds(20);
  uint32_t v = 0, cnt = 1<<TOUCH_Y_N;
  do{v+=adc_single_read(ADC_TOUCH_Y);}while(--cnt);
  return v>>TOUCH_Y_N;
}

static int
touch_measure_x(void)
{
  // drive high to low on Y line (coordinates from left to right)
  palSetPad(GPIOB, GPIOB_YN);
  palClearPad(GPIOA, GPIOA_YP);
  // Set Y line as output
  palSetPadMode(GPIOB, GPIOB_YN, PAL_MODE_OUTPUT_PUSHPULL);
  palSetPadMode(GPIOA, GPIOA_YP, PAL_MODE_OUTPUT_PUSHPULL);
  // Set X line as input
  palSetPadMode(GPIOB, GPIOB_XN, PAL_MODE_INPUT);        // Hi-z mode
  palSetPadMode(GPIOA, GPIOA_XP, PAL_MODE_INPUT_ANALOG); // <- ADC_TOUCH_X channel

  uint32_t v = 0, cnt = 1<<TOUCH_X_N;
  do{v+=adc_single_read(ADC_TOUCH_X);}while(--cnt);
  return v>>TOUCH_X_N;
}
// Manually measure touch event
static inline int
touch_status(void)
{
  return adc_single_read(ADC_TOUCH_Y) > TOUCH_THRESHOLD;
}

static void
touch_prepare_sense(void)
{
  // Set Y line as input
  palSetPadMode(GPIOB, GPIOB_YN, PAL_MODE_INPUT);          // Hi-z mode
  palSetPadMode(GPIOA, GPIOA_YP, PAL_MODE_INPUT_PULLDOWN); // Use pull
  // drive high on X line (for touch sense on Y)
  palSetPad(GPIOB, GPIOB_XN);
  palSetPad(GPIOA, GPIOA_XP);
  // force high X line
  palSetPadMode(GPIOB, GPIOB_XN, PAL_MODE_OUTPUT_PUSHPULL);
  palSetPadMode(GPIOA, GPIOA_XP, PAL_MODE_OUTPUT_PUSHPULL);

//  chThdSleepMilliseconds(10); // Wait 10ms for denounce touch
}

static void
touch_start_watchdog(void)
{
  if (touch_status_flag&TOUCH_INTERRUPT_ENABLED) return;
  touch_status_flag^=TOUCH_INTERRUPT_ENABLED;
  adc_start_analog_watchdog();
#ifdef __REMOTE_DESKTOP__
  remote_mouse_down = 0;
#endif
}
static void
touch_stop_watchdog(void)
{
  if (!(touch_status_flag&TOUCH_INTERRUPT_ENABLED)) return;
  touch_status_flag^=TOUCH_INTERRUPT_ENABLED;
  adc_stop_analog_watchdog();
}

// Touch panel timer check (check press frequency 20Hz)
static const GPTConfig gpt3cfg = {
  20,     // 200Hz timer clock. 200/10 = 20Hz touch check
  NULL,   // Timer callback.
  0x0020, // CR2:MMS=02 to output TRGO
  0
};

//
// Touch init function init timer 3 trigger adc for check touch interrupt, and run measure
//
static void touch_init(void){
  // Prepare pin for measure touch event
  touch_prepare_sense();
  // Start touch interrupt, used timer_3 ADC check threshold:
  gptStart(&GPTD3, &gpt3cfg);         // Init timer 3
  gptStartContinuous(&GPTD3, 10);     // Start timer 3 vs timer 10 interval
  touch_start_watchdog();             // Start ADC watchdog (measure by timer 3 interval and trigger interrupt if touch pressed)
}

// Main software touch function, should:
// set last_touch_x and last_touch_x
// return touch status
static int
touch_check(void)
{
  touch_stop_watchdog();

  int stat = touch_status();
  if (stat) {
    int y = touch_measure_y();
    int x = touch_measure_x();
    touch_prepare_sense();
    if (touch_status())
    {
      last_touch_x = x;
      last_touch_y = y;
    }
#ifdef __REMOTE_DESKTOP__
    remote_mouse_down = 0;
  } else {
    stat = remote_mouse_down == 1;
#endif
  }
  if (stat != last_touch_status) {
    last_touch_status = stat;
    return stat ? EVT_TOUCH_PRESSED : EVT_TOUCH_RELEASED;
  }
  return stat ? EVT_TOUCH_DOWN : EVT_TOUCH_NONE;
}
//*******************************************************************************
// End Software Touch module
//*******************************************************************************
#endif // end SOFTWARE_TOUCH

void touch_set(int16_t x, int16_t y) {
  last_touch_x = x;
  last_touch_y = y;
}

void
touch_wait_release(void)
{
  while (touch_check() != EVT_TOUCH_NONE)
    chThdSleepMilliseconds(20);
}
#if 0
static inline void
touch_wait_pressed(void)
{
  while (touch_check() != EVT_TOUCH_PRESSED)
    ;
}
#endif

static inline void
touch_wait_released(void)
{
  while (touch_check() != EVT_TOUCH_RELEASED)
    ;
}

void
touch_cal_exec(void)
{
  int x1, x2, y1, y2;
  int old_flip = config.flip;
  config.flip = 0;
  ili9341_set_foreground(LCD_FG_COLOR);
  ili9341_set_background(LCD_BG_COLOR);
  ili9341_clear_screen();
  ili9341_line(0, 0, 0, 32);
  ili9341_line(0, 0, 32, 0);
  ili9341_line(0, 0, 32, 32);
  ili9341_drawstring("TOUCH UPPER LEFT", 40, 40);
  touch_wait_released();
//  touch_wait_release();
  x1 = last_touch_x;
  y1 = last_touch_y;

  ili9341_clear_screen();
  ili9341_line(LCD_WIDTH-1, LCD_HEIGHT-1, LCD_WIDTH-1, LCD_HEIGHT-32);
  ili9341_line(LCD_WIDTH-1, LCD_HEIGHT-1, LCD_WIDTH-32, LCD_HEIGHT-1);
  ili9341_line(LCD_WIDTH-1, LCD_HEIGHT-1, LCD_WIDTH-32, LCD_HEIGHT-32);
  ili9341_drawstring("TOUCH LOWER RIGHT", LCD_WIDTH-17*(FONT_WIDTH)-30, LCD_HEIGHT-FONT_GET_HEIGHT-35);

  touch_wait_released();
 // touch_wait_release();
  x2 = last_touch_x;
  y2 = last_touch_y;

  config.touch_cal[0] = x1;
  config.touch_cal[1] = y1;
  config.touch_cal[2] = (x2 - x1) * 16 / LCD_WIDTH;
  config.touch_cal[3] = (y2 - y1) * 16 / LCD_HEIGHT;

  config.flip = old_flip;

  config_save();            // Auto save touch calibration

  //redraw_all();
}

void
touch_draw_test(void)
{
  int x0, y0;
  int x1, y1;

  ili9341_set_foreground(LCD_FG_COLOR);
  ili9341_set_background(LCD_BG_COLOR);
  ili9341_clear_screen();
  ili9341_drawstring("TOUCH TEST: DRAG PANEL, PRESS BUTTON TO FINISH", OFFSETX, LCD_HEIGHT - FONT_GET_HEIGHT);

  int old_button_state = 0;
  while (touch_check() != EVT_TOUCH_PRESSED) {
    int button_state = READ_PORT() & BUTTON_MASK;
    if (button_state != old_button_state) {
      char buf[20];
      plot_printf(buf, sizeof buf, "STATE: %4d       ", button_state);
      ili9341_drawstring_7x13(buf, 120, 120);
      old_button_state = button_state;
    }

  }

  do {
    if (touch_check() == EVT_TOUCH_PRESSED){
      touch_position(&x0, &y0);
      do {
        chThdSleepMilliseconds(50);
        touch_position(&x1, &y1);
        ili9341_line(x0, y0, x1, y1);
        x0 = x1;
        y0 = y1;
      } while (touch_check() != EVT_TOUCH_RELEASED);
    }
  }while (!(btn_check() & EVT_BUTTON_SINGLE_CLICK));
}


void
touch_position(int *x, int *y)
{
#ifdef __REMOTE_DESKTOP__
  if (remote_mouse_down) {
    *x = last_touch_x;
    *y = last_touch_y;
    return;
  }
#endif
  int tx = (last_touch_x - config.touch_cal[0]) * 16 / config.touch_cal[2];
  int ty = (last_touch_y - config.touch_cal[1]) * 16 / config.touch_cal[3];

  if (config.flip) {
    tx = LCD_WIDTH - 1 - tx;
    ty = LCD_HEIGHT - 1 - ty;
  }
  *x = tx;
  *y = ty;
}
#ifdef TINYSA4
extern const char *get_hw_version_text(void);
#endif

void
show_version(void)
{
  int x = 5, y = 5, i = 0;
  ili9341_set_foreground(LCD_FG_COLOR);
  ili9341_set_background(LCD_BG_COLOR);

  ili9341_clear_screen();
  uint16_t shift = 0b00000100001;
// Version text for tinySA3
#ifdef TINYSA3
  ili9341_drawstring_10x14(info_about[i++], x , y);
  y+=FONT_GET_HEIGHT*3+3-5;
  while (info_about[i]) {
    do {shift>>=1; y+=5;} while (shift&1);
    ili9341_drawstring(info_about[i++], x, y+=FONT_STR_HEIGHT+3-5);
  }
  if (has_esd)
    ili9341_drawstring("ESD protected", x, y+=FONT_STR_HEIGHT + 2);

  y+=FONT_STR_HEIGHT + 1;
#endif
// Version text for tinySA4
#ifdef TINYSA4
  ili9341_set_foreground(LCD_MENU_ACTIVE_COLOR);
  ili9341_drawstring_10x14(info_about[i++], x , y);
  ili9341_set_foreground(LCD_FG_COLOR);
  y+=FONT_GET_HEIGHT*3+2-5;
  ili9341_drawstring_7x13(info_about[i++], x , y);
  while (info_about[i]) {
    do {shift>>=1; y+=5;} while (shift&1);
    ili9341_drawstring_7x13(info_about[i++], x, y+=bFONT_STR_HEIGHT+2-5);
  }
  char buf[96];
  plot_printf(buf, sizeof(buf), "HW Version:%s", get_hw_version_text());
  ili9341_drawstring_7x13(buf, x, y+=bFONT_STR_HEIGHT);

extern const char *states[];
#define ENABLE_THREADS_COMMAND
#ifdef ENABLE_THREADS_COMMAND
  y+=FONT_STR_HEIGHT + 1;
  thread_t *tp;
  tp = chRegFirstThread();
  do {
    uint32_t max_stack_use = 0U;
#if (CH_DBG_ENABLE_STACK_CHECK == TRUE) || (CH_CFG_USE_DYNAMIC == TRUE)
    uint32_t stklimit = (uint32_t)tp->wabase;
#if CH_DBG_FILL_THREADS == TRUE
    uint8_t *p = (uint8_t *)tp->wabase; while(p[max_stack_use]==CH_DBG_STACK_FILL_VALUE) max_stack_use++;
#endif
#else
    uint32_t stklimit = 0U;
#endif
    char buf[96];
    plot_printf(buf, sizeof(buf), "%08x|%08x|%08x|%08x|%4u|%4u|%9s|%12s",
             stklimit, (uint32_t)tp->ctx.sp, max_stack_use, (uint32_t)tp,
             (uint32_t)tp->refs - 1, (uint32_t)tp->prio, states[tp->state],
             tp->name == NULL ? "" : tp->name);
    ili9341_drawstring_7x13(buf, x, y+=bFONT_STR_HEIGHT);
    tp = chRegNextThread(tp);
  } while (tp != NULL);
#endif
  y+=bFONT_STR_HEIGHT + 1;
#endif  // TINYSA4
  uint16_t cnt = 0;
  while (true) {
    if (touch_check() == EVT_TOUCH_PRESSED)
      break;
    if (btn_check() & EVT_BUTTON_SINGLE_CLICK)
      break;
    chThdSleepMilliseconds(40);
    if ((cnt++)&0x07) continue; // Not update time so fast

#ifdef TINYSA4
#ifdef __USE_RTC__
    uint32_t tr = rtc_get_tr_bin(); // TR read first
    uint32_t dr = rtc_get_dr_bin(); // DR read second
    char buf[96];
    plot_printf(buf, sizeof(buf), "Time: 20%02d/%02d/%02d %02d:%02d:%02d" " (LS%c)",
      RTC_DR_YEAR(dr),
      RTC_DR_MONTH(dr),
      RTC_DR_DAY(dr),
      RTC_TR_HOUR(dr),
      RTC_TR_MIN(dr),
      RTC_TR_SEC(dr),
      (RCC->BDCR & STM32_RTCSEL_MASK) == STM32_RTCSEL_LSE ? 'E' : 'I');
    ili9341_drawstring_7x13(buf, x, y);
#endif
#if 0
    uint32_t vbat=adc_vbat_read();
    plot_printf(buf, sizeof(buf), "Batt: %d.%03dV", vbat/1000, vbat%1000);
    ili9341_drawstring_7x13(buf, x, y + bFONT_STR_HEIGHT + 1);
#endif
#endif // TINYSA4
  }
}

#ifndef TINYSA4
void
enter_dfu(void)
{
  int x = 5, y = 5;
  ili9341_set_foreground(LCD_FG_COLOR);
  ili9341_set_background(LCD_BG_COLOR);
  // leave a last message 
  ili9341_clear_screen();
  ili9341_drawstring_7x13("DFU: Device Firmware Update Mode\n"
                          "To exit DFU mode, please reset device yourself.", x, y);
  // see __early_init in ./NANOVNA_STM32_F072/board.c
  *((unsigned long *)BOOT_FROM_SYTEM_MEMORY_MAGIC_ADDRESS) = BOOT_FROM_SYTEM_MEMORY_MAGIC;
  NVIC_SystemReset();
}
#endif

static void
select_lever_mode(int mode)
{
  if (uistat.lever_mode != mode) {
    uistat.lever_mode = mode;
    redraw_request |= REDRAW_FREQUENCY | REDRAW_MARKER;
  }
}

// type of menu item 
enum {
  MT_NONE,                      // sentinel menu
//  MT_BLANK,                     // blank menu (nothing draw)
  MT_SUBMENU,                   // enter to submenu
  MT_CALLBACK,                  // call user function
  MT_ADV_CALLBACK,              // adv call user function
  MT_CANCEL,                    // menu, step back on one level up
  MT_TITLE,                     // Title
  MT_KEYPAD,
  MT_REPEATS = 0x08,
  MT_ICON = 0x10,
  MT_HIGH = 0x20,               // Only applicable to high mode
  MT_LOW = 0x40,                // Only applicable to low mode
  MT_FORM = 0x80,               // Large button menu
};
//#define MT_BACK     0x40
//#define MT_LEAVE    0x20
#define MT_MASK(x) (0x7 & (x))
#define DATA_STARTS_REPEATS(S,R)  (((R)<<4)+(S))

#define MT_CUSTOM_LABEL  0

// Call back functions for MT_CALLBACK type
typedef void (*menuaction_cb_t)(int item, uint16_t data);
#define UI_FUNCTION_CALLBACK(ui_function_name) void ui_function_name(int item, uint16_t data)

typedef void (*menuaction_acb_t)(int item, uint16_t data, ui_button_t *b);
#define UI_FUNCTION_ADV_CALLBACK(ui_function_name) void ui_function_name(int item, uint16_t data, ui_button_t *b)

static freq_t
get_marker_frequency(int marker)
{
  if (marker < 0 || marker >= MARKERS_MAX)
    return 0;
  if (!markers[marker].enabled)
    return 0;
  return getFrequency(markers[marker].index);
}

static UI_FUNCTION_CALLBACK(menu_marker_op_cb)
{
  (void)item;
  freq_t freq = get_marker_frequency(active_marker);
  if (freq == 0)
    return; // no active marker

  switch (data) {
  case 0: /* MARKER->START */
  case 1: /* MARKER->STOP */
  case 2: /* MARKER->CENTER */
    set_sweep_frequency(data, freq);
    if (data == 2) {
      uistat.lever_mode = LM_SPAN;
      uistat.auto_center_marker = true;
    }
    break;
  case 3: /* MARKERS->SPAN */
    {
      if (previous_marker == MARKER_INVALID || active_marker == previous_marker) {
        // if only 1 marker is active, keep center freq and make span the marker comes to the edge
        freq_t center = get_sweep_frequency(ST_CENTER);
        freq_t span = center > freq ? center - freq : freq - center;
        set_sweep_frequency(ST_SPAN, span * 2);
      } else {
        // if 2 or more marker active, set start and stop freq to each marker
        freq_t freq2 = get_marker_frequency(previous_marker);
        if (freq2 == 0)
          return;
        if (freq > freq2) {
          freq2 = freq;
          freq = get_marker_frequency(previous_marker);
        }
        set_sweep_frequency(ST_START, freq);
        set_sweep_frequency(ST_STOP, freq2);
      }
    }
    break;
  case 4: // marker -> ref level
    {
    float l = actual_t[markers[active_marker].index];
    float s_max = value(l)/setting.scale;
    user_set_reflevel(setting.scale*(floorf(s_max)+2));
    }
    break;
#ifdef __VNA__
  case 4: /* MARKERS->EDELAY */
    { 
      if (uistat.current_trace == -1)
        break;
      float (*array)[2] = measured[trace[uistat.current_trace].channel];
      float v = groupdelay_from_array(markers[active_marker].index, array);
      set_electrical_delay(electrical_delay + (v / 1e-12));
    }
    break;
#endif
  }
  menu_move_back(true);
  redraw_request |= REDRAW_CAL_STATUS;
  //redraw_all();
}

static UI_FUNCTION_CALLBACK(menu_markers_reset_cb)
{
  (void)item;
  (void)data;
  markers_reset();
}

static UI_FUNCTION_CALLBACK(menu_marker_search_cb)
{
  (void)item;
  int i = -1;
  if (active_marker == MARKER_INVALID)
    return;
  markers[active_marker].mtype &= ~M_TRACKING;
  switch (data) {
  case 0: /* search Left */
    i = marker_search_left_min(active_marker);
    break;
  case 1: /* search right */
    i = marker_search_right_min(active_marker);
    break;
  case 2: /* search Left */
    i = marker_search_left_max(active_marker);
    break;
  case 3: /* search right */
    i = marker_search_right_max(active_marker);
    break;
  case 4: /* peak search */
    i = marker_search_max(active_marker);
    break;
  }
  if (i != -1) {
    markers[active_marker].index = i;
    if (data > 1) // Maximum related
      interpolate_maximum(active_marker);
    else
      markers[active_marker].frequency = getFrequency(i);
  }
  redraw_marker(active_marker);
//  if (data == 4)
    select_lever_mode(LM_MARKER);   // Allow any position with level
//  else
//    select_lever_mode(LM_SEARCH); // Jump from maximum to maximum
}

#if 0
static UI_FUNCTION_ADV_CALLBACK(menu_marker_tracking_acb){
  (void)item;
  (void)data;
  if (active_marker == MARKER_INVALID) return;
  if(b){
    b->icon = markers[active_marker].mtype & M_TRACKING ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  markers[active_marker].mtype ^= M_TRACKING;
}
#endif

#ifdef __VNA__
static void
menu_marker_smith_cb(int item, uint8_t data)
{
  (void)item;
  marker_smith_format = data;
  redraw_marker(active_marker);
  draw_menu();
}
#endif

static void
active_marker_select(int item)  // used only to select an active marker from the modify marker selection menu
{
  if (item == -1) {
    active_marker = previous_marker;
    previous_marker = MARKER_INVALID;
    if (active_marker == MARKER_INVALID) {
      choose_active_marker();
    }
  } else {
    if (previous_marker != active_marker) {
      previous_marker = active_marker;
      active_marker = item;
    } else {
      active_marker = item;
    }
  }
}

//----------------------- iu_sa.c inserted -------------------

#define FORM_ICON_WIDTH      16
#define FORM_ICON_HEIGHT     16

static const uint8_t left_icons [] =
{
#define I_EMPTY   (0*16)
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),

#define I_HIGH_INPUT (1*16)
                              // +----------------+
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000001100000), // |         **     |
  _BMP16(0b0000000000111001), // |          ***  *|
  _BMP16(0b0000111111111111), // |   *************|
  _BMP16(0b0000000000111001), // |          ***  *|
  _BMP16(0b0000000001100000), // |         **     |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000001), // |                |
  _BMP16(0b0000000000000001), // |                |
  _BMP16(0b0000000000000001), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
                              // +----------------+
#define I_LOW_INPUT (2*16)
                              // +----------------+
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000001), // |                |
  _BMP16(0b0000000000000001), // |                |
  _BMP16(0b0000000000000001), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000001100000), // |         **     |
  _BMP16(0b0000000000111001), // |         ****  *|
  _BMP16(0b0000111111111111), // |   *************|
  _BMP16(0b0000000000111001), // |         ****  *|
  _BMP16(0b0000000001100000), // |         **     |
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0000000000000000), // |                |
                              // +----------------+
#define I_LOW_OUTPUT (3*16)
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000110000000),
  _BMP16(0b0000011100000001),
  _BMP16(0b0000111111111111),
  _BMP16(0b0000011100000001),
  _BMP16(0b0000000110000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),

#define I_HIGH_OUTPUT (4*16)
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000110000000),
  _BMP16(0b0000011100000001),
  _BMP16(0b0000111111111111),
  _BMP16(0b0000011100000001),
  _BMP16(0b0000000110000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000001),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),

#define I_CONNECT (5*16)
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000110000),
  _BMP16(0b0000000000111101),
  _BMP16(0b0000001111111111),
  _BMP16(0b0000010000111101),
  _BMP16(0b0000100000110000),
  _BMP16(0b0001000000000000),
  _BMP16(0b0001000000000000),
  _BMP16(0b0000100000110000),
  _BMP16(0b0000010000111101),
  _BMP16(0b0000001111111111),
  _BMP16(0b0000000000111101),
  _BMP16(0b0000000000110000),
  _BMP16(0b0000000000000000),
  _BMP16(0b0000000000000000),
};

const uint8_t right_icons [] =
{
#define I_SA    0
                              // +----------------+
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0111111111111111), // | ***************|
  _BMP16(0b0100000000000001), // | *             *|
  _BMP16(0b1100000000000001), // |**             *|
  _BMP16(0b1100000000000001), // |**  *          *|
  _BMP16(0b1100000000000001), // |**  *          *|
  _BMP16(0b0100100000000001), // | *  *   *      *|
  _BMP16(0b0100100000000001), // | *  *   *      *|
  _BMP16(0b0100101010001001), // | *  *   *   *  *|
  _BMP16(0b0100101010101001), // | *  * * *   *  *|
  _BMP16(0b1100101010101001), // |**  * * * * *  *|
  _BMP16(0b1101111111111101), // |**  * * * * *  *|
  _BMP16(0b1100000000000001), // |** *********** *|
  _BMP16(0b0100000000000001), // | *             *|
  _BMP16(0b0111111111111111), // | ***************|
  _BMP16(0b0000000000000000), // |                |
                              // +----------------+

#define I_GEN   1
                              // +----------------+
  _BMP16(0b0000000000000000), // |                |
  _BMP16(0b0111111111111111), // | ***************|
  _BMP16(0b0100000000000001), // | *             *|
  _BMP16(0b1100000000000001), // |**             *|
  _BMP16(0b1100111110001101), // |**  *****   ** *|
  _BMP16(0b1100100010001001), // |**  *   *   *  *|
  _BMP16(0b0100100010001001), // | *  *   *   *  *|
  _BMP16(0b0100100010001001), // | *  *   *   *  *|
  _BMP16(0b0100100010001001), // | *  *   *   *  *|
  _BMP16(0b0100100010001001), // | *  *   *   *  *|
  _BMP16(0b1100100010001001), // |**  *   *   *  *|
  _BMP16(0b1101100011111001), // |** **   *****  *|
  _BMP16(0b1100000000000001), // |**             *|
  _BMP16(0b0100000000000001), // | *             *|
  _BMP16(0b0111111111111111), // | ***************|
  _BMP16(0b0000000000000000), // |                |
                              // +----------------+
#define I_CONFIG 2
  _BMP16(0b0000000000000000),
  _BMP16(0b0111111111111111),
  _BMP16(0b0100000000000001),
  _BMP16(0b1100000010000001),
  _BMP16(0b1100001111000001),
  _BMP16(0b1100011110001001),
  _BMP16(0b0100011100011101),
  _BMP16(0b0100011110111001),
  _BMP16(0b0100001111111001),
  _BMP16(0b0100011111110001),
  _BMP16(0b1100111110000001),
  _BMP16(0b1101111100000001),
  _BMP16(0b1100111000000001),
  _BMP16(0b0100000000000001),
  _BMP16(0b0111111111111111),
  _BMP16(0b0000000000000000),

#define I_SINUS 3
  _BMP16(0b0000000000000000),
  _BMP16(0b0111111111111111),  // 1
  _BMP16(0b0100000000000001),  // 2
  _BMP16(0b1100000000000001),  // 3
  _BMP16(0b1100000000110001),  // 4
  _BMP16(0b1100000001001001),  // 5
  _BMP16(0b0100000010000101),  // 6
  _BMP16(0b0101000010000101),  // 7
  _BMP16(0b0101000010000101),  // 8
  _BMP16(0b0101000010000001),  // 9
  _BMP16(0b1100100100000001),  //10
  _BMP16(0b1100011000000001),  //11
  _BMP16(0b1100000000000001),  //12
  _BMP16(0b0100000000000001),  //13
  _BMP16(0b0111111111111111),  //14
  _BMP16(0b0000000000000000),
};

#define KP_X(x) (48*(x) + 2 + (LCD_WIDTH-BUTTON_WIDTH-192))
#define KP_Y(y) (48*(y) + 2)


#define KP_PERIOD 10
#define KP_MINUS 11
#define KP_X1 12
#define KP_K 13
#define KP_M 14
#define KP_G 15
#define KP_BS 16
#define KP_INF 17
#define KP_DB 18
#define KP_PLUSMINUS 19
#define KP_KEYPAD 20
#define KP_m 21
#define KP_u 22
#define KP_n 23
#define KP_p 24

#define KP_0    31
#define KP_1    32
#define KP_2    33
#define KP_5    34
#define KP_10   35
#define KP_20   36
#define KP_50   37
#define KP_100  38
#define KP_200  39
#define KP_500  40


typedef struct {
  uint8_t x:4;
  uint8_t y:4;
  int8_t  c;
} keypads_t;

static const keypads_t *keypads;

// 7 8 9 G
// 4 5 6 M
// 1 2 3 k
// 0 . < x

static const keypads_t keypads_freq[] = {
  { 1, 3, KP_PERIOD },
  { 0, 3, 0 },
  { 0, 2, 1 },
  { 1, 2, 2 },
  { 2, 2, 3 },
  { 0, 1, 4 },
  { 1, 1, 5 },
  { 2, 1, 6 },
  { 0, 0, 7 },
  { 1, 0, 8 },
  { 2, 0, 9 },
  { 3, 0, KP_G },
  { 3, 1, KP_M },
  { 3, 2, KP_K },
  { 3, 3, KP_X1 },
  { 2, 3, KP_BS },
  { 0, 0, -1 }
};

// 7 8 9
// 4 5 6
// 1 2 3
// 0 . < x

static const keypads_t keypads_positive[] = {
  { 1, 3, KP_PERIOD },
  { 0, 3, 0 },
  { 0, 2, 1 },
  { 1, 2, 2 },
  { 2, 2, 3 },
  { 0, 1, 4 },
  { 1, 1, 5 },
  { 2, 1, 6 },
  { 0, 0, 7 },
  { 1, 0, 8 },
  { 2, 0, 9 },
  { 3, 3, KP_X1 },
  { 2, 3, KP_BS },
  { 0, 0, -1 }
};

// 100 200 500 n
// 10  20  50  u
// 1   2   5   m
// 0   .   <   x

static const keypads_t keypads_pos_unit[] = {
  { 1, 3, KP_PERIOD },
  { 0, 3, 0 },
  { 0, 2, 1 },
  { 1, 2, 2 },
  { 2, 2, 5 },
  { 0, 1, KP_10 },
  { 1, 1, KP_20 },
  { 2, 1, KP_50 },
  { 0, 0, KP_100 },
  { 1, 0, KP_200 },
  { 2, 0, KP_500 },
  { 3, 0, KP_n },
  { 3, 1, KP_u },
  { 3, 2, KP_m },
  { 3, 3, KP_X1 },
  { 2, 3, KP_BS },
  { 0, 0, -1 }
};

// 7 8 9 m
// 4 5 6 u
// 1 2 3 -
// 0 . < x

static const keypads_t keypads_plusmin_unit[] = {
  { 1, 3, KP_PERIOD },
  { 0, 3, 0 },
  { 0, 2, 1 },
  { 1, 2, 2 },
  { 2, 2, 3 },
  { 0, 1, 4 },
  { 1, 1, 5 },
  { 2, 1, 6 },
  { 0, 0, 7 },
  { 1, 0, 8 },
  { 2, 0, 9 },
  { 3, 0, KP_u},
  { 3, 1, KP_m},
  { 3, 2, KP_MINUS },
  { 3, 3, KP_X1 },
  { 2, 3, KP_BS },
  { 0, 0, -1 }
};
// 7 8 9
// 4 5 6
// 1 2 3 -
// 0 . < x

static const keypads_t keypads_plusmin[] = {
  { 1, 3, KP_PERIOD },
  { 0, 3, 0 },
  { 0, 2, 1 },
  { 1, 2, 2 },
  { 2, 2, 3 },
  { 0, 1, 4 },
  { 1, 1, 5 },
  { 2, 1, 6 },
  { 0, 0, 7 },
  { 1, 0, 8 },
  { 2, 0, 9 },
  { 3, 0, KP_u},
  { 3, 1, KP_m},
  { 3, 2, KP_MINUS },
  { 3, 3, KP_X1 },
  { 2, 3, KP_BS },
  { 0, 0, -1 }
};

// 7 8 9
// 4 5 6
// 1 2 3 m
// 0 . < x
static const keypads_t keypads_time[] = {
  { 1, 3, KP_PERIOD },
  { 0, 3, 0 },
  { 0, 2, 1 },
  { 1, 2, 2 },
  { 2, 2, 3 },
  { 0, 1, 4 },
  { 1, 1, 5 },
  { 2, 1, 6 },
  { 0, 0, 7 },
  { 1, 0, 8 },
  { 2, 0, 9 },
//  { 3, 0, KP_n},
//  { 3, 1, KP_u},
  { 3, 2, KP_m },
  { 3, 3, KP_X1 },
  { 2, 3, KP_BS },
  { 0, 0, -1 }
};

enum {
  KM_START, KM_STOP, KM_CENTER, KM_SPAN, KM_CW, // These must be first to share common help text
  //#5
  KM_REFLEVEL, KM_SCALE, KM_ATTENUATION, KM_ACTUALPOWER, KM_IF,
  // #10
  KM_SAMPLETIME, KM_LOWOUTLEVEL, KM_DECAY, KM_NOISE,
#ifdef TINYSA4
  KM_FREQ_CORR,
#else
  KM_10MHZ,
#endif
  // #15
  KM_REPEAT, KM_EXT_GAIN, KM_TRIGGER, KM_LEVELSWEEP, KM_SWEEP_TIME,
  // #20
  KM_OFFSET_DELAY, KM_FAST_SPEEDUP, KM_GRIDLINES, KM_MARKER, KM_MODULATION,
  // #25
  KM_HIGHOUTLEVEL,
#ifdef TINYSA4
  KM_COR_AM,  KM_COR_WFM, KM_COR_NFM,
  KM_DEVIATION, KM_DEPTH,
  KM_IF2,
  // #30
  KM_R,KM_MOD,KM_CP,
#endif
  KM_ATTACK,
#ifdef __ULTRA__
  KM_ULTRA_START,
#endif
#ifdef TINYSA4
  KM_EXP_AVER,
#endif
  KM_LEVEL,
#ifdef __LIMITS__
  KM_LIMIT_FREQ, KM_LIMIT_LEVEL,
#endif
  KM_MARKER_TIME,
  // #35
  KM_VAR,
#ifdef __NOISE_FIGURE__
  KM_NF,
#endif
  KM_LINEAR_SCALE,
#ifdef TINYSA4
  KM_DIRECT_START,
  KM_DIRECT_STOP,
#ifdef __USE_RTC__
  KM_RTC_DATE,
  KM_RTC_TIME,
#endif
#endif
  KM_CODE,
  KM_NONE // always at enum end
};

static const struct {
  const keypads_t *keypad_type;
  char * name;
} keypads_mode_tbl[KM_NONE] = {
[KM_START]        = {keypads_freq        , "START"}, // start
[KM_STOP]         = {keypads_freq        , "STOP"}, // stop
[KM_CENTER]       = {keypads_freq        , "CENTER"}, // center
[KM_SPAN]         = {keypads_freq        , "SPAN"}, // span
[KM_CW]           = {keypads_freq        , "FREQ"}, // cw freq
[KM_REFLEVEL]     = {keypads_plusmin_unit, "REF\nLEVEL"}, // reflevel #5
[KM_SCALE]        = {keypads_pos_unit    , "SCALE"}, // scale
[KM_ATTENUATION]  = {keypads_positive    , "ATTENUATE"}, // attenuation
[KM_ACTUALPOWER]  = {keypads_plusmin_unit, "ACTUAL\nPOWER"}, // actual power
[KM_IF]           = {keypads_freq        , "IF"}, // IF
[KM_SAMPLETIME]   = {keypads_positive    , "SAMPLE\nDELAY"}, // sample delay #10
[KM_LOWOUTLEVEL]  = {keypads_plusmin     , "OUTPUT\nLEVEL"},    // KM_LOWOUTLEVEL
[KM_DECAY]        = {keypads_positive    , "DECAY"},    // KM_DECAY
[KM_NOISE]        = {keypads_positive    , "NOISE\nLEVEL"},    // KM_NOISE
#ifdef TINYSA4
[KM_FREQ_CORR]    = {keypads_plusmin     , "PPB"},    // KM_FREQ_CORR
#else
[KM_10MHZ]        = {keypads_freq        , "FREQ"},    // KM_10MHz
#endif
[KM_REPEAT]       = {keypads_positive    , "SAMPLE\nREPEAT"},    // KM_REPEA #15
[KM_EXT_GAIN]     = {keypads_plusmin     , "EXT\nGAIN"},    // KM_EXT_GAIN
[KM_TRIGGER]      = {keypads_plusmin_unit, "LEVEL"},    // KM_TRIGGER
[KM_LEVELSWEEP]   = {keypads_plusmin     , "LEVEL\nSWEEP"},    // KM_LEVELSWEEP
[KM_SWEEP_TIME]   = {keypads_time        , "SWEEP\nSECONDS"},    // KM_SWEEP_TIME
[KM_OFFSET_DELAY] = {keypads_positive    , "OFFSET\nDELAY"}, // KM_OFFSET_DELAY #20
[KM_FAST_SPEEDUP] = {keypads_positive    , "FAST\nSPEEDUP"}, // KM_FAST_SPEEDUP
[KM_GRIDLINES]    = {keypads_positive    , "MINIMUM\nGRIDLINES"}, // KM_GRIDLINES
[KM_MARKER]       = {keypads_freq        , "MARKER\nFREQ"}, // KM_MARKER
[KM_MODULATION]   = {keypads_freq        , "MODULATION\nFREQ"}, // KM_MODULATION
[KM_HIGHOUTLEVEL] = {keypads_plusmin     , "OUTPUT\nLEVEL"},    // KM_HIGHOUTLEVEL #25
#ifdef TINYSA4
[KM_COR_AM]       = {keypads_plusmin     , "COR\nAM"},    // KM_COR_AM
[KM_COR_WFM]      = {keypads_plusmin     , "COR\nWFM"},    // KM_COR_WFM
[KM_COR_NFM]      = {keypads_plusmin     , "COR\nNFM"},    // KM_COR_NFM
[KM_DEVIATION]    = {keypads_freq        , "DEVIATION"}, // KM_DEVIATION
[KM_DEPTH]        = {keypads_positive    , "DEPTH%"}, // KM_DEPTH
[KM_IF2]          = {keypads_freq        , "IF2"}, // KM_IF2
[KM_R]            = {keypads_plusmin     , "R"}, // KM_R    #30
[KM_MOD]          = {keypads_positive    , "MODULO"}, // KM_MOD
[KM_CP]           = {keypads_positive    , "CP"}, // KM_CP
#endif
[KM_ATTACK]       = {keypads_positive    , "ATTACK"},    // KM_ATTACK
#ifdef __ULTRA__
[KM_ULTRA_START]  = {keypads_freq        , "ULTRA\nSTART"}, // KM_ULTRA_START
#endif
#ifdef TINYSA4
[KM_EXP_AVER]     = {keypads_positive    , "EXPONENTIAL\nAVERAGING"}, //KM_EXP_AVER
#endif
[KM_LEVEL]        = {keypads_plusmin     , "LEVEL"}, // KM_LEVEL
#ifdef __LIMITS__
[KM_LIMIT_FREQ]   = {keypads_freq         , "FREQ"},  // KM_LIMIT_FREQ
[KM_LIMIT_LEVEL]  = {keypads_plusmin_unit , "LEVEL"},  // KM_LIMIT_LEVEL
#endif
[KM_MARKER_TIME]  = {keypads_time        , "MARKER\nTIME"}, // KM_MARKER_TIME
[KM_VAR]          = {keypads_freq        , "JOG\nSTEP"}, // jog step
#ifdef __NOISE_FIGURE__
[KM_NF]           = {keypads_plusmin     , "NOISE\nFIGURE"}, // noise figure of tinySA
#endif
[KM_LINEAR_SCALE] = {keypads_plusmin     , "SCALE"}, // scale for linear units
#ifdef TINYSA4
[KM_DIRECT_START] = {keypads_freq        , "DIRECT\nSTART"}, // KM_DIRECT_START
[KM_DIRECT_STOP]  = {keypads_freq        , "DIRECT\nSTOP"}, // KM_DIRECT_STOP
#ifdef __USE_RTC__
[KM_RTC_DATE]     = {keypads_positive    , "SET DATE\n YYMMDD"}, // Date
[KM_RTC_TIME]     = {keypads_positive    , "SET TIME\n HHMMSS"}, // Time
#endif
#endif
[KM_CODE]         = {keypads_positive    , "CODE"},              // KM_CODE
};

#if 0 // Not used
enum { SL_GENLOW_FREQ, SL_GENHIGH_FREQ, SL_GENLOW_LEVEL, SL_GENHIGH_LEVEL };
ui_slider_t ui_sliders [] =
{
 { KM_CENTER,       true, 0, 1000000,   0,          350000000,  M_GENLOW},
 { KM_CENTER,       true, 0, 1000000,   240000000,  960000000,  M_GENHIGH},
 { KM_LOWOUTLEVEL,  false,0, 1,         -76,        -6,         M_GENLOW},
 { KM_HIGHOUTLEVEL, false,0, 1,         -38,        +6,         M_GENHIGH},
};
#endif


// ===[MENU CALLBACKS]=========================================================
const menuitem_t  menu_lowoutputmode[];
const menuitem_t  menu_highoutputmode[];
const menuitem_t  menu_mode[];
static const menuitem_t  menu_modulation[];
static const menuitem_t  menu_top[];
static const menuitem_t  menu_trace[];
static const menuitem_t  menu_marker_trace[];
static const menuitem_t  menu_subtract_trace[];
#ifdef __LIMITS__
static const menuitem_t  menu_limit_modify[];
static const menuitem_t  menu_limit_select[];
#endif
static const menuitem_t  menu_average[];
static const menuitem_t  menu_reffer[];
static const menuitem_t  menu_sweep_points[];
static const menuitem_t  menu_sweep_points_form[];
static const menuitem_t  menu_modulation[];
static const menuitem_t  menu_marker_ref_select[];
#ifdef __USE_SERIAL_CONSOLE__
static const menuitem_t  menu_connection[];
#endif
//static const menuitem_t  menu_drive_wide[];
static const menuitem_t  menu_config[];
#ifdef TINYSA4
static const menuitem_t  menu_settings3[];
static const menuitem_t  menu_curve[];
static const menuitem_t  menu_curve_confirm[];
static const menuitem_t  menu_measure_noise_figure[];
static const menuitem_t  menu_calibrate_harmonic[];
#endif
static const menuitem_t  menu_sweep[];
static const menuitem_t  menu_settings[];
static const menuitem_t  menu_settings2[];
static const menuitem_t  menu_lowoutput_settings[];
extern bool dirty;
char range_text[20];

#ifdef TINYSA4
int input_is_calibrated(void)
{
  if (config.input_is_calibrated)
    return true;
  drawMessageBox("Error", "First calibrate 100kHz to 5.34GHz input", 2000);
  redraw_request|= REDRAW_AREA;
  return false;
}

int output_is_calibrated(void)
{
  if (config.output_is_calibrated)
    return true;
  drawMessageBox("Error", "First calibrate 30MHz output", 2000);
  redraw_request|= REDRAW_AREA;
  return false;
}

static int unlock_internals = 0;

static UI_FUNCTION_ADV_CALLBACK(menu_internals_acb)
{
  (void)data;
  (void)item;
  if (b){
    return;
  }
  if (unlock_internals != 5432) {
    kp_help_text = "Internals access code";
    ui_mode_keypad(KM_CODE);
    if (uistat.value != 5432) {
      return;
    }
    unlock_internals = 5432;
  }
  menu_push_submenu(menu_settings2);
}

#endif

static UI_FUNCTION_ADV_CALLBACK(menu_sweep_acb)
{
  (void)data;
  (void)item;
  if (b){
    if (setting.level_sweep != 0 || get_sweep_frequency(ST_SPAN) != 0)  {
      plot_printf(b->text, sizeof b->text, "SW:%3.2fMHz %+ddB %.3Fs",
                  get_sweep_frequency(ST_SPAN) / 1000000.0,
                  (int)setting.level_sweep,
                  setting.sweep_time_us/(float)ONE_SECOND_TIME);
    }
    else
      plot_printf(b->text, sizeof b->text, "SWEEP: OFF");
    return;
  }
  menu_push_submenu(menu_sweep);
}

#ifdef __SWEEP_RESTART__
static UI_FUNCTION_ADV_CALLBACK(menu_restart_acb){
  (void)item;
  (void)data;
  if(b){
    if (setting.sweep) {
      if (current_index >= 0) {
        float current_level = setting.level + ((float)current_index)* setting.level_sweep / (float)sweep_points;
        plot_printf(b->text, sizeof b->text, "STOP %5.3QHz %+.1fdBm", getFrequency(current_index), current_level);
      } else
        plot_printf(b->text, sizeof b->text, "STOP SWEEP");
    }
    else
      plot_printf(b->text, sizeof b->text, "START SWEEP");
    return;
  }
  setting.sweep = !setting.sweep;
  dirty = true;
}
#endif

#ifdef TINYSA4

float local_actual_level;
int current_curve;
int current_curve_index;

static UI_FUNCTION_ADV_CALLBACK(menu_curve_acb)
        {
  (void)item;
  int old_m;
  if (b){
    plot_printf(b->text, sizeof b->text, "%8.3QHz %+4.1fdB",
                config.correction_frequency[current_curve][data],
                config.correction_value[current_curve][data]);
    return;
  }
  switch(current_curve) {
  case CORRECTION_LOW_OUT:
    old_m = setting.mode;
    reset_settings(M_GENLOW);
    set_level(-35);
    set_sweep_frequency(ST_CW, config.correction_frequency[current_curve][data]);
    setting.mute = false;
    perform(false, 0, config.correction_frequency[current_curve][data], false);
    perform(false, 1, config.correction_frequency[current_curve][data], false);
    plot_printf(uistat.text, sizeof uistat.text, "Level of %.3QHz output",
                config.correction_frequency[current_curve][data]);
    kp_help_text = uistat.text;
    kp_buf[0]=0;
    ui_mode_keypad(KM_LEVEL);

    if (kp_buf[0] != 0) {
      float new_offset = (-35.0) - uistat.value + config.correction_value[current_curve][data];        // calculate offset based on difference between measured peak level and known peak level
      if (new_offset > -25 && new_offset < 25) {
        config.correction_value[current_curve][data] = new_offset;
        config_save();
      }
    }
    reset_settings(old_m);
    break;
  case CORRECTION_LNA:
    reset_settings(M_LOW);
    setting.extra_lna = true;
    goto common;
  case CORRECTION_LOW:
    reset_settings(M_LOW);
    common:
    set_sweep_frequency(ST_SPAN,   1000000);
    set_sweep_frequency(ST_CENTER, config.correction_frequency[current_curve][data]);
    setting.step_delay_mode = SD_PRECISE;
    current_curve_index = data;
    menu_push_submenu(menu_curve_confirm);
    break;
  }
}

extern float peakLevel;

UI_FUNCTION_CALLBACK(menu_curve_confirm_cb)
{
  (void)item;
  if (data) {
    float new_offset = local_actual_level - peakLevel + config.correction_value[current_curve][current_curve_index];        // calculate offset based on difference between measured peak level and known peak level
    if (new_offset > -30 && new_offset < 30) {
      config.correction_value[current_curve][current_curve_index] = new_offset;
      config_save();
    }
  }
  menu_move_back(false);
}

float measured_noise_figure;

#if 0
UI_FUNCTION_CALLBACK(menu_noise_figure_confirm_cb)
{
  (void)item;
  if (data) {
    if (measured_noise_figure > 3 && measured_noise_figure < 15) {
      config.noise_figure = measured_noise_figure;
      config_save();
      nf_gain = 0.00001;                            // almost zero
      set_measurement(M_NF_VALIDATE);               // Continue to validate
      return;
    }
  }
  menu_move_back(false);
}
#endif

static UI_FUNCTION_CALLBACK(menu_input_curve_prepare_cb)
{
  (void)item;
  (void)data;
  if (!input_is_calibrated())
    return;
  kp_help_text = "Enter actual input level";
  kp_buf[0]=0;
  ui_mode_keypad(KM_LEVEL);
  if (kp_buf[0] != 0) {
    local_actual_level = uistat.value;
    current_curve = CORRECTION_LOW;
    menu_push_submenu(menu_curve);
  }
}

static UI_FUNCTION_CALLBACK(menu_lna_curve_prepare_cb)
{
  (void)item;
  (void)data;
  if (!input_is_calibrated())
    return;
  kp_help_text = "Enter actual input level";
  kp_buf[0]=0;
  ui_mode_keypad(KM_LEVEL);
  if (kp_buf[0] != 0) {
    local_actual_level = uistat.value;
    current_curve = CORRECTION_LNA;
    menu_push_submenu(menu_curve);
  }
}

static UI_FUNCTION_CALLBACK(menu_lna_u_curve_prepare_cb)
{
  (void)item;
  (void)data;
  if (!input_is_calibrated())
    return;
  kp_help_text = "Enter actual input level";
  kp_buf[0]=0;
  ui_mode_keypad(KM_LEVEL);
  if (kp_buf[0] != 0) {
    local_actual_level = uistat.value;
    current_curve = CORRECTION_LNA_ULTRA;
    menu_push_submenu(menu_curve);
  }
}

static UI_FUNCTION_CALLBACK(menu_ultra_curve_prepare_cb)
{
  (void)item;
  (void)data;
  if (!input_is_calibrated())
    return;
  kp_help_text = "Enter actual input level";
  kp_buf[0]=0;
  ui_mode_keypad(KM_LEVEL);
  if (kp_buf[0] != 0) {
    local_actual_level = uistat.value;
    current_curve = CORRECTION_LOW_ULTRA;
    menu_push_submenu(menu_curve);
  }
}

static UI_FUNCTION_CALLBACK(menu_output_curve_prepare_cb)
{
  (void)item;
  (void)data;
  if (!output_is_calibrated())
    return;
  current_curve = CORRECTION_LOW_OUT;
  menu_push_submenu(menu_curve);
}

#endif

static UI_FUNCTION_ADV_CALLBACK(menu_output_level_acb)
{
  (void)item;
  (void)data;
  if (b){
    return;
  }
  int old_m = setting.mode;
  reset_settings(M_GENLOW);
#ifdef TINYSA4
#define TEST_LEVEL -30
#else
#define TEST_LEVEL -25
#endif
  set_level(TEST_LEVEL);
  set_sweep_frequency(ST_CW, 30000000);
  setting.mute = false;
  perform(false, 0, 30000000, false);
  perform(false, 1, 30000000, false);
  kp_help_text = "Enter actual level of 30MHz output";
  kp_buf[0]=0;
  ui_mode_keypad(KM_LEVEL);
  if (kp_buf[0] != 0) {
    float old_offset = config.low_level_output_offset;
    if (!config.input_is_calibrated) old_offset = 0;
    float new_offset = uistat.value - (TEST_LEVEL) + old_offset;        // calculate offset based on difference between measured peak level and known peak level
    if (uistat.value == 100) { new_offset = 0; config.input_is_calibrated = false; }
    if (new_offset > -15 && new_offset < 15) {
      config.output_is_calibrated = true;
      config.low_level_output_offset = new_offset;
      config_save();
    }
  }
  reset_settings(old_m);
}

#ifdef TINYSA4
static UI_FUNCTION_ADV_CALLBACK(menu_output_level2_acb)
{
  (void)item;
  (void)data;
  if (b){
    return;
  }
  if (!output_is_calibrated())
    return;
  int old_m = setting.mode;
  reset_settings(M_GENLOW);
  set_level(-30);
  set_sweep_frequency(ST_CW, 1000000000);
  setting.mute = false;
  perform(false, 0, 1000000000, false);
  perform(false, 1, 1000000000, false);
  kp_help_text = "Enter actual level of 1GHz output";
  kp_buf[0]=0;
  ui_mode_keypad(KM_LEVEL);
  if (kp_buf[0] != 0) {
    float old_offset = config.direct_level_output_offset;
    float new_offset = (-30.0) - uistat.value  + old_offset;        // calculate offset based on difference between measured peak level and known peak level
    if (new_offset > -10 && new_offset < 10) {
      config.direct_level_output_offset = new_offset;
      config_save();
    }
  }
  reset_settings(old_m);
}

static UI_FUNCTION_ADV_CALLBACK(menu_output_level3_acb)
{
  (void)item;
  (void)data;
  if (b){
    return;
  }
  if (!output_is_calibrated())
    return;
  int old_m = setting.mode;
  reset_settings(M_GENLOW);

  force_signal_path = true;
  test_path = PATH_LEAKAGE;
  test_output_drive = -1;

  set_level(-30);
  set_sweep_frequency(ST_CW, 1200000000);
  setting.mute = false;
  perform(false, 0, 1200000000, false);
  perform(false, 1, 1200000000, false);
  kp_help_text = "Enter actual level of 1.2GHz output";
  kp_buf[0]=0;
  ui_mode_keypad(KM_LEVEL);
  if (kp_buf[0] != 0) {
    float old_offset = config.adf_level_offset;
    float new_offset = (-30.0) - uistat.value  + old_offset;        // calculate offset based on difference between measured peak level and known peak level
    if (new_offset > -10 && new_offset < 10) {
      config.adf_level_offset = new_offset;
      config_save();
    }
  }
  force_signal_path = false;
  reset_settings(old_m);
}

#endif

#ifdef TINYSA4
static const int item_to_mode[2] = { 0,2 };
#else
static const int item_to_mode[4] = { 0,1,2,3 };
#endif
static UI_FUNCTION_ADV_CALLBACK(menu_mode_acb)
{
  (void)data;
  item = item_to_mode[item];
  if (b){
    if (item == setting.mode)  {
      b->param_1.text = "Return";
      b->bg = LCD_MENU_ACTIVE_COLOR;
      b->border = BUTTON_BORDER_FALLING | MENU_BUTTON_BORDER;
    }
    else
      b->param_1.text = "Switch";
    return;
  }
  set_mode(item);
//  draw_cal_status();
  switch (item) {
  case 0:
//    if (setting.mode != M_LOW)
//      set_mode(M_LOW);
    menu_move_back(true);
    break;
  case 1:
//    if (setting.mode != M_HIGH)
//      set_mode(M_HIGH);
    menu_move_back(true);
    break;
  case 2:
    menu_push_submenu(menu_lowoutputmode);
    break;
  case 3:
    menu_push_submenu(menu_highoutputmode);
    break;
  }
  redraw_request |= REDRAW_CAL_STATUS;
}

static UI_FUNCTION_ADV_CALLBACK(menu_load_preset_acb)
{
  (void)item;
  if(b){
    setting_t *p = caldata_pointer(data);
    if (p)
      plot_printf(b->text, sizeof(b->text), "%.6FHz\n%.6FHz", (float)p->frequency0, (float)p->frequency1);
    else
      plot_printf(b->text, sizeof(b->text), "EMPTY %d", (int)data);
    return;
  }
  if (caldata_recall(data) == -1) {
    if (data == 0)
      reset_settings(setting.mode);  // Restore factory defaults
  }
  menu_move_back(true);
}

static UI_FUNCTION_ADV_CALLBACK(menu_store_preset_acb)
{
  (void)item;
  if(b){
    b->param_1.u = data;
    return;
  }
  if (data == 100) {
    reset_settings(M_LOW);  // Restore all defaults in Low mode
    set_refer_output(-1);
 //   setting.mode = -1;
    data = 0;
  }
  caldata_save(data);
  menu_move_back(true);
}

#ifdef __SD_CARD_LOAD__
UI_FUNCTION_CALLBACK(menu_load_config_cb)
{
  (void)item;
  (void)data;
  sd_card_load_config("config.ini");
  ui_mode_normal();
}
#endif

UI_FUNCTION_CALLBACK(menu_autosettings_cb)
{
  (void)item;
  (void)data;
  reset_settings(setting.mode);

  markers_reset();
  //  set_refer_output(1);

  //  SetPowerLevel(100); // Reset
//  set_clear_storage();
  dirty = true;
  //  menu_move_back(true);   // stay in input menu
  ui_mode_normal();
//  draw_cal_status();
}

static UI_FUNCTION_CALLBACK(menu_scale_cb)
{
  (void)item;
  (void)data;
  kp_help_text = "Enter scale";
  kp_buf[0]=0;

  if (UNIT_IS_LINEAR(setting.unit))
    ui_mode_keypad(KM_LINEAR_SCALE);
  else
    ui_mode_keypad(KM_SCALE);
  ui_mode_normal();
}

#ifdef __CALIBRATE__
static UI_FUNCTION_CALLBACK(menu_calibrate_cb)
{
  (void)data;
  (void)item;
  switch (data) {
  case 1:
    sweep_mode = SWEEP_CALIBRATE;
#ifdef TINYSA4
    menu_move_back(false);
#endif
    menu_move_back(true);
    break;
  case 2:
    reset_calibration();
    break;
#ifdef TINYSA4
  case 3:
    if (!input_is_calibrated())
      return;
    sweep_mode = SWEEP_CALIBRATE_HARMONIC;
    menu_move_back(false);
    menu_move_back(true);
    break;
#endif
  }
}
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_scanning_speed_acb)
{
  (void)item;
  if(b){
    b->icon = data == setting.step_delay_mode ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  set_step_delay(data);
//    menu_move_back(false);
  ui_mode_normal();
}

#define CONFIG_MENUITEM_TOUCH_CAL   0
#define CONFIG_MENUITEM_TOUCH_TEST  1
#define CONFIG_MENUITEM_SELFTEST    2
#define CONFIG_MENUITEM_VERSION     3
static UI_FUNCTION_CALLBACK(menu_config_cb)
{
  (void)item;
  switch (data) {
  case CONFIG_MENUITEM_TOUCH_CAL:
    touch_cal_exec();
    break;
  case CONFIG_MENUITEM_TOUCH_TEST:
    touch_draw_test();
    break;
  case CONFIG_MENUITEM_SELFTEST:
    sweep_mode = 0;         // Suspend sweep to save time
    menu_move_back(true);
    setting.test = 0;
    setting.test_argument = 0;
    sweep_mode = SWEEP_SELFTEST;
    return;
  case CONFIG_MENUITEM_VERSION:
    show_version();
    break;
  }
  ui_mode_normal();
  redraw_frame();
  request_to_redraw_grid();
}
#ifndef TINYSA4
static UI_FUNCTION_CALLBACK(menu_dfu_cb)
{
  (void)data;
  (void)item;
  enter_dfu();
}
#endif

#ifdef __LISTEN__
static UI_FUNCTION_ADV_CALLBACK(menu_listen_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = (sweep_mode & SWEEP_LISTEN) ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  if (sweep_mode & SWEEP_LISTEN) {
    sweep_mode = SWEEP_ENABLE;
  } else {
    sweep_mode = SWEEP_LISTEN;
  }
  ui_mode_normal();
  redraw_frame();
  request_to_redraw_grid();

#if 0
  if (markers[active_marker].enabled == M_ENABLED) {
    do {
      perform(false,0,frequencies[markers[active_marker].index], false);
      SI4432_Listen(MODE_SELECT(setting.mode));
    } while (ui_process_listen_lever());
  }
#endif
}
#endif
#ifdef TINYSA4


static UI_FUNCTION_ADV_CALLBACK(menu_lowoutput_settings_acb)
{
  static char mode_string[26];
  (void)item;
  if (b){
    if (data == 255) {
      plot_printf(mode_string, sizeof mode_string, "%s %s %s %s",
                  (!force_signal_path ? "" : path_text[test_path]),
                  (get_sweep_frequency(ST_START) < MINIMUM_DIRECT_FREQ ? "SINUS" : "" ),
                  (get_sweep_frequency(ST_STOP) >= MINIMUM_DIRECT_FREQ ? "SQUARE WAVE" : ""),
                  (get_sweep_frequency(ST_STOP) > MAX_LOW_OUTPUT_FREQ && setting.mixer_output ? "MIXER" : ""));
      b->param_1.text = mode_string;
      return; }
    b->icon = data == setting.mixer_output ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  switch(data) {
  case 255:
    menu_push_submenu(menu_lowoutput_settings);
    return;
  case 0:
    setting.mixer_output = false;
    dirty = true;
    break;
  case 1:
    setting.mixer_output = true;
    dirty = true;
    break;
  }
  menu_move_back(false);
}

#endif
// const int menu_modulation_value[]={MO_NONE,MO_AM, MO_NFM, MO_WFM, MO_EXTERNAL};
const char *menu_modulation_text[MO_MAX]=
{  "None", "AM 30%",
#ifdef TINYSA4
   "FM 2.5kHz",
   "FM 3kHz",
   "FM 5kHz",
#else
   "FM 4kHz",
#endif
   "FM 75kHz", "External"};

static UI_FUNCTION_ADV_CALLBACK(menu_modulation_acb)
{
  (void)item;
  if (b){
    plot_printf(b->text, sizeof b->text, "%s", menu_modulation_text[data]);
    b->icon = data == setting.modulation ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
//Serial.println(item);
  if (data) {
    set_sweep_frequency(ST_SPAN, 0);      // No other scanning allowed when modulation is on!!!!!
    set_level_sweep(0);
  }
  set_modulation(data);
//  menu_move_back(false);  // Don't move back
}

static UI_FUNCTION_ADV_CALLBACK(menu_smodulation_acb){
  (void)item;
  (void)data;
  if(b){
    if (setting.modulation == MO_NONE || setting.modulation == MO_EXTERNAL)
      plot_printf(b->text, sizeof b->text, "MOD: %s", menu_modulation_text[setting.modulation]);
    else {
#ifdef TINYSA4
      if (setting.modulation == MO_AM)
        plot_printf(b->text, sizeof b->text, "MOD: %4dHz AM %d%%", (int)(setting.modulation_frequency), setting.modulation_depth_x100);
      else
        plot_printf(b->text, sizeof b->text, "MOD: %4dHz FM %4QHz", (int)(setting.modulation_frequency), (freq_t)(setting.modulation_deviation_div100*100));
#else
      plot_printf(b->text, sizeof b->text, "MOD: %4dHz %s", (int)(setting.modulation_frequency), menu_modulation_text[setting.modulation]);
#endif
    }
    return;
  }
  menu_push_submenu(menu_modulation);
}

//                               0      1       2       3      4      5      6      7
const char *menu_reffer_text[]={"OFF","30MHz","15MHz","10MHz","4MHz","3MHz","2MHz","1MHz"};
static UI_FUNCTION_ADV_CALLBACK(menu_reffer_acb)
{
  (void)item;
  if (b){
    b->icon = setting.refer == ((int)data-1) ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    b->param_1.text = menu_reffer_text[data];
    return;
  }
//Serial.println(item);
  set_refer_output((int)data - 1);
  menu_move_back(false);
//  ui_mode_normal();   // Stay in menu mode
//  draw_cal_status();
}

static UI_FUNCTION_ADV_CALLBACK(menu_sreffer_acb){
  (void)item;
  (void)data;
  if(b){
    b->param_1.text = menu_reffer_text[setting.refer+1];
    return;
  }
  menu_push_submenu(menu_reffer);
}

#ifdef TINYSA3
static UI_FUNCTION_ADV_CALLBACK(menu_lo_drive_acb)
{
  (void)item;
  if(b){
    b->param_1.i = drive_dBm[data];
    b->icon = data == setting.lo_drive ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
//Serial.println(item);
  set_lo_drive(data);
  menu_move_back(false);
//  ui_mode_normal();
//  draw_cal_status();
}
#else
static UI_FUNCTION_ADV_CALLBACK(menu_mixer_drive_acb)
{
  (void)item;
  if(b){
    b->param_1.i = data;
    b->icon = data == setting.lo_drive ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
//Serial.println(item);
  set_lo_drive(data);
  menu_move_back(false);
//  ui_mode_normal();
//  draw_cal_status();
}
#endif

#if 0
static UI_FUNCTION_ADV_CALLBACK(menu_sdrive_acb){
  (void)item;
  (void)data;
  if(b){
#ifdef TINYSA4
    b->param_1.i = setting.lo_drive;
#else
    b->param_1.i = drive_dBm[setting.lo_drive] + (setting.mode==M_GENHIGH ? setting.external_gain : 0);
#endif
    return;
  }
  menu_push_submenu(menu_drive_wide);
}
#endif

#ifdef __SPUR__
static UI_FUNCTION_ADV_CALLBACK(menu_spur_acb)
{
  (void)data;
  (void)item;
  if (b){
    if (setting.mode == M_LOW) {
      b->param_1.text = "SPUR\nREMOVAL";
      b->icon = AUTO_ICON(setting.spur_removal);
    } else {
      b->param_1.text = "MIRROR\nMASKING";
#ifdef TINYSA4
      b->icon = AUTO_ICON(setting.mirror_masking ? 1 : 0);  // mirror_masking does not yet have an auto mode so this is never an auto icon
#else
      b->icon = setting.mirror_masking == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
#endif
    }
    return;
  }
  if (setting.mode == M_LOW) {
    toggle_spur();
  } else
    toggle_mirror_masking();
  //  menu_move_back(false);
  ui_mode_normal();
}
#if 0
#ifdef __HARMONIC__
static UI_FUNCTION_ADV_CALLBACK(menu_harmonic_spur_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = AUTO_ICON(setting.spur_removal);
    return;
  }
  toggle_spur();
  ui_mode_normal();
}#endif
#endif
#endif
#endif

#ifdef __ULTRA__
static UI_FUNCTION_ADV_CALLBACK(menu_debug_spur_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = debug_spur == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  toggle_debug_spur();
  //  menu_move_back();
  ui_mode_normal();
}
#endif

#ifdef TINYSA4
static UI_FUNCTION_ADV_CALLBACK(menu_hide_21MHz_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = config.hide_21MHz == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  config.hide_21MHz = ! config.hide_21MHz;
  //  menu_move_back();
  ui_mode_normal();
}

static UI_FUNCTION_ADV_CALLBACK(menu_progress_bar_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = progress_bar == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  progress_bar = !progress_bar;
  //  menu_move_back();
  ui_mode_normal();
}

static UI_FUNCTION_ADV_CALLBACK(menu_extra_lna_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = setting.extra_lna == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  toggle_extra_lna();
  //  menu_move_back(false);
  ui_mode_normal();
}
#ifndef __NEW_SWITCHES__
static UI_FUNCTION_ADV_CALLBACK(menu_adf_out_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = config.high_out_adf4350 == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  toggle_high_out_adf4350();
  //  menu_move_back(false);
  ui_mode_normal();
}
#endif

#ifdef __WAIT_CTS_WHILE_SLEEPING__
volatile int sleep = 0;
static UI_FUNCTION_ADV_CALLBACK(menu_sleep_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = sleep == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  sleep = !sleep;
}
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_debug_avoid_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = debug_avoid == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  toggle_debug_avoid();
  //  menu_move_back();
  ui_mode_normal();
}

static UI_FUNCTION_ADV_CALLBACK(menu_debug_freq_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = debug_frequencies == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  debug_frequencies = ! debug_frequencies;
  //  menu_move_back();
  ui_mode_normal();
}


static UI_FUNCTION_ADV_CALLBACK(menu_linear_averaging_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = linear_averaging == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  linear_averaging = ! linear_averaging;
  dirty = true;
  //  menu_move_back();
  ui_mode_normal();
}


#endif

#ifdef __ULTRA__
static UI_FUNCTION_ADV_CALLBACK(menu_ultra_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = config.ultra == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  if (!config.ultra) {
    drawMessageBox("Info", "Visit tinysa.org/ultra for unlock code", 2000);

    kp_help_text = "Ultra unlock code";
    ui_mode_keypad(KM_CODE);
    if (uistat.value != 4321)
      return;
  }
  config.ultra = !config.ultra;
  config_save();
  reset_settings(M_LOW);
  if (config.ultra){
    set_sweep_frequency(ST_START, 0);
    set_sweep_frequency(ST_STOP, 3000000000ULL);
  }
  //  menu_move_back(false);
  ui_mode_normal();
}

static UI_FUNCTION_ADV_CALLBACK(menu_direct_acb)
{
  (void)data;
  (void)item;
  if (b){
    b->icon = config.direct == 0 ? BUTTON_ICON_NOCHECK : BUTTON_ICON_CHECK;
    return;
  }
  config.direct = !config.direct;
  config_save();
  //  menu_move_back(false);
  ui_mode_normal();
}

#endif

static UI_FUNCTION_CALLBACK(menu_clearconfig_cb)
{
  (void)data;
  (void)item;
  kp_help_text = "Clear unlock code";
  ui_mode_keypad(KM_CODE);
  if (uistat.value != 1234)
    return;
  clear_all_config_prop_data();
#ifndef TINYSA4
#define SCB_AIRCR_VECTKEYSTAT_LSB 16
#define SCB_AIRCR_VECTKEYSTAT       (0xFFFF << SCB_AIRCR_VECTKEYSTAT_LSB)
#define SCB_AIRCR_VECTKEY       (0x05FA << SCB_AIRCR_VECTKEYSTAT_LSB)
#define SCB_AIRCR_SYSRESETREQ         (1 << 2)
  SCB->AIRCR = SCB_AIRCR_VECTKEY | SCB_AIRCR_SYSRESETREQ;
          while(true);
#else
  rccEnableWWDG(FALSE);
  WWDG->CFR = 0x60;
  WWDG->CR = 0xff;
  /* wait forever */
#endif
  while (1)
    ;
}

float nf_gain;
const char * const averageText[] = { "OFF", "MINH", "MAXH", "MAXD", " AVER4", "A16", "AVER", "QUASI", "TABLE", "DECONV"};

static UI_FUNCTION_ADV_CALLBACK(menu_measure_acb)
{
  (void)item;
  if (b){
    b->icon = data == setting.measurement ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  menu_move_back(false);
  markers_reset();

  if ((data != M_OFF && setting.measurement != M_OFF) || data == M_OFF )
  {
    //      reset_settings(setting.mode);
    if (0) {
    no_measurement:
      drawMessageBox("Error", "Incorrect input", 2000);
      redraw_request|= REDRAW_AREA;

      data = M_OFF;
    }
    if (setting.measurement == M_LINEARITY) {
      TRACE_DISABLE(TRACE_STORED_FLAG);
    }
    set_average(0,AV_OFF);
    set_external_gain(0.0);
#ifdef TINYSA4
    set_extra_lna(false);
#endif
  }
  switch(data) {
    case M_OFF:                                     // Off
//      set_measurement(M_OFF);
      break;
    case M_IMD:                                     // IMD
      reset_settings(setting.mode);
      for (int i = 1; i< MARKERS_MAX; i++) {
        markers[i].mtype = M_DELTA;
      }
      kp_help_text = "Frequency of fundamental";
      ui_mode_keypad(KM_CENTER);
      set_sweep_frequency(ST_START, 0);
      set_sweep_frequency(ST_STOP, uistat.value*(MARKERS_MAX+1));
      set_average(0,AV_4);
//      set_measurement(M_IMD);
      break;
    case M_OIP3:                                     // OIP3
      reset_settings(setting.mode);
      for (int i = 0; i< 4; i++) {
        markers[i].enabled = M_ENABLED;
      }
      markers[1].mtype = M_TRACKING;
      kp_help_text = "Frequency of left signal";
      ui_mode_keypad(KM_CENTER);
      freq_t left =  uistat.value;
      kp_help_text = "Right signal";
      ui_mode_keypad(KM_CENTER);
      freq_t right =  uistat.value;
      set_sweep_frequency(ST_CENTER, (left+right)/2);
      freq_t local_span = (right - left)*10;
#ifdef TINYSA4
      if (local_span < 3000)
        local_span = 3000;
#else
      if (local_span < 30000)
        local_span = 30000;
#endif
      set_sweep_frequency(ST_SPAN, local_span);
      set_average(0,AV_4);
//      set_measurement(M_OIP3);
      break;
    case M_PHASE_NOISE:                             // Phase noise
      reset_settings(setting.mode);
      markers[1].enabled = M_ENABLED;
      markers[1].mtype = M_DELTA | M_NOISE;
      kp_help_text = "Frequency of signal";
      ui_mode_keypad(KM_CENTER);
      kp_help_text = "Frequency offset";
      ui_mode_keypad(KM_SPAN);
      set_sweep_frequency(ST_SPAN, uistat.value*4);
//      set_measurement(M_PHASE_NOISE);
      set_average(0,AV_4);

      break;
    case M_SNR:                             // STop band measurement
      reset_settings(setting.mode);
      kp_help_text = "Frequency of signal";
      ui_mode_keypad(KM_CENTER);
      kp_help_text = "Bandwidth";
      ui_mode_keypad(KM_SPAN);
      set_sweep_frequency(ST_SPAN, uistat.value*3);
//      set_measurement(M_SNR);
      set_average(0,AV_4);

      break;
    case M_PASS_BAND:                             // pass band measurement
//      reset_settings(setting.mode);
      markers[1].enabled = M_ENABLED;
      markers[2].enabled = M_ENABLED;
//      kp_help_text = "Frequency of signal";
//      ui_mode_keypad(KM_CENTER);
//      kp_help_text = "Width of signal";
//      ui_mode_keypad(KM_SPAN);
//      set_sweep_frequency(ST_SPAN, uistat.value*2);
      set_measurement(M_PASS_BAND);
//      SetAverage(4);

      break;
#ifdef __LINEARITY__
    case M_LINEARITY:
      TRACE_ENABLE(TRACE_STORED_FLAG);
//      set_measurement(M_LINEARITY);
      break;
#endif
    case M_AM:                                     // AM
      reset_settings(setting.mode);
      for (int i = 1; i< 3; i++) {
        markers[i].enabled = M_ENABLED;
#ifdef TINYSA4
        markers[i].mtype = M_DELTA| M_TRACKING;
#else
//        markers[i].mtype = M_DELTA;// | M_TRACKING;
#endif
      }
#ifdef TINYSA4
      freq_t span;
#else
      freq_t center, span;
      markers[0].mtype = M_NORMAL; // M_REFERENCE | M_TRACKING;           // Not M_TRACKING!!!!
#endif
      kp_help_text = "Frequency of signal";
      ui_mode_keypad(KM_CENTER);
#ifdef TINYSA3
      center = uistat.value;
#endif
#ifdef TINYSA4
      kp_help_text = "Modulation frequency: 500hz .. 10kHz";
#else
      kp_help_text = "Modulation frequency: 3 .. 10kHz";
#endif
      ui_mode_keypad(KM_SPAN);
      span = uistat.value;
#ifdef TINYSA4
      if (span < 500 || span > 10000)
        goto no_measurement;
#endif
#ifdef TINYSA4
      set_RBW((span * 5 / 50) / 100);
#endif
      set_sweep_frequency(ST_SPAN, span * 5);
//      update_frequencies();                     // To ensure markers are positioned right!!!!!!
//      set_measurement(M_AM);
#ifndef TINYSA4
      set_marker_frequency(0, center);
      set_marker_frequency(1, center-span);
      set_marker_frequency(2, center+span);
#endif
      set_average(0,AV_4);
      break;
    case M_FM:                                     // FM
      reset_settings(setting.mode);
      for (int i = 1; i< 3; i++) {
        markers[i].enabled = M_ENABLED;
      }
      markers[0].mtype = M_NORMAL;              /// Not M_TRACKING !!!!
      kp_help_text = "Frequency of signal";
      ui_mode_keypad(KM_CENTER);
      set_marker_frequency(0, uistat.value);
#ifdef TINYSA4
      kp_help_text = "Modulation frequency: 500Hz .. 10kHz";
      ui_mode_keypad(KM_SPAN);
      if (uistat.value < 500 || uistat.value > 10000)
        goto no_measurement;
      set_RBW(uistat.value/300);
#else
      kp_help_text = "Modulation frequency: 1 .. 2.5kHz";
      ui_mode_keypad(KM_SPAN);
      if (uistat.value < 1000 || uistat.value > 2500)
        goto no_measurement;
      set_RBW(uistat.value/100);
#endif
      // actual_rbw_x10
#ifdef TINYSA4
      kp_help_text = "Frequency deviation: 3 .. 500kHz";
#define MINIMUM_DEVIATION   1500
#else
      kp_help_text = "Frequency deviation: 500Hz .. 500kHz";
#define MINIMUM_DEVIATION   12000
#endif
      ui_mode_keypad(KM_SPAN);
      if (uistat.value < MINIMUM_DEVIATION)
        uistat.value = MINIMUM_DEVIATION;   // minimum span
      set_sweep_frequency(ST_SPAN, uistat.value*4);
//      set_measurement(M_FM);
      break;
    case M_THD:
      set_average(0,AV_4);
//      set_measurement(M_THD);
      break;
#ifdef __CHANNEL_POWER__
    case M_CP:                             // channel power
      reset_settings(setting.mode);
      markers[0].enabled = M_DISABLED;
      kp_help_text = "Channel frequency";
      ui_mode_keypad(KM_CENTER);
      kp_help_text = "Channel width";
      ui_mode_keypad(KM_SPAN);
      set_sweep_frequency(ST_SPAN, uistat.value*3);
//      set_measurement(M_CP);
      break;
#endif
#ifdef __NOISE_FIGURE__
    case M_NF_TINYSA:
      reset_settings(setting.mode);
      set_refer_output(-1);
      nf_gain = 0;
      goto noise_figure;
    case M_NF_STORE:
      if (measured_noise_figure > 2 && measured_noise_figure < 20) {
        config.noise_figure = measured_noise_figure;
        config_save();
        data = M_NF_VALIDATE;               // Continue to validate
        goto validate;
      } else {
        drawMessageBox("Error", "NF out of range",1000);
        data = M_NF_TINYSA;               // Continue to measure
      }
      break;
    case M_NF_VALIDATE:
validate:
      nf_gain = 0.00001;                            // almost zero
      goto noise_figure;
    case M_NF_AMPLIFIER:                             // noise figure
#if 0
      reset_settings(setting.mode);
      set_refer_output(-1);
#endif
      kp_help_text = "Amplifier Gain ";
      float old_gain = setting.external_gain;
      ui_mode_keypad(KM_EXT_GAIN);
      nf_gain = setting.external_gain;
      setting.external_gain = old_gain;
  noise_figure:
      markers[0].enabled = M_ENABLED;
      markers[0].mtype = M_NOISE | M_AVER;          // Not tracking
      set_extra_lna(true);
      set_attenuation(0);
      if (data == M_NF_TINYSA) {
        kp_help_text = "Noise center frequency";
        ui_mode_keypad(KM_CENTER);
        set_marker_frequency(0, uistat.value);
#if 0
        kp_help_text = "Noise span";
        ui_mode_keypad(KM_SPAN);
#else
        set_sweep_frequency(ST_SPAN, 100000);
#endif
        set_RBW(1000);  // 300kHz
      }

//      set_sweep_frequency(ST_SPAN, 0);
      set_average(0,AV_100);
      if (data == M_NF_TINYSA || data == M_NF_VALIDATE ) {
        menu_push_submenu(menu_measure_noise_figure);
        goto leave;
      }
      break;
#endif
#ifdef __FFT_DECONV__
    case M_DECONV:
      set_average(0,AV_DECONV);
      break;
#endif
  }

//  selection = -1;
  ui_mode_normal();
  goto leave;       // to get rid of warning
leave:
  set_measurement(data);
//  draw_cal_status();
}

static UI_FUNCTION_ADV_CALLBACK(menu_atten_acb)
{
  (void)item;
  (void)data;
  if(b){
    if  (MODE_HIGH(setting.mode)) {
      b->param_1.text = "0dB";
      b->icon = (setting.atten_step*30 == data) ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    }
    else {
      b->param_1.text = "AUTO";
      b->icon = setting.auto_attenuation ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    }
    return;
  }
  if  (MODE_HIGH(setting.mode)) {
    setting.auto_attenuation = false;
    set_attenuation(0);
  } else
    set_auto_attenuation();
  menu_move_back(true);
}

static UI_FUNCTION_ADV_CALLBACK(menu_atten_high_acb)
{
  (void)item;
  if(b){
    b->icon = (setting.atten_step*30 == data) ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  setting.auto_attenuation = false;
  set_attenuation(data);
  menu_move_back(true);
}

static UI_FUNCTION_ADV_CALLBACK(menu_reflevel_acb)
{
  (void)item;
  (void)data;
  if(b){
    b->icon = setting.auto_reflevel ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  set_auto_reflevel(true);
  menu_move_back(true);
}

static uint8_t current_trace = TRACE_ACTUAL;

static UI_FUNCTION_ADV_CALLBACK(menu_average_acb)
{
  (void)item;
  if (b){
    b->icon = setting.average[current_trace] == data ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  set_average(current_trace,data);
  if (data == AV_TABLE)
    menu_push_submenu(menu_limit_select);
  else {
    ui_mode_normal();
  }
//  menu_move_back(true);
}

static UI_FUNCTION_ADV_CALLBACK(menu_trace_acb)
{
  (void)item;
  if(b){
    b->param_1.i = data+1;
    b->icon = (data == current_trace) ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    b->bg = LCD_TRACE_1_COLOR+data;
    return;
  }
  if (setting.normalized_trace != -1 && data == TRACE_TEMP) {
    drawMessageBox("Error", "Disable normalization first", 2000);
    redraw_request|= REDRAW_AREA;
  } else
    current_trace = data;
  menu_move_back(false);
}

static UI_FUNCTION_ADV_CALLBACK(menu_marker_trace_acb)
{
  (void)item;
  if(b){
    b->param_1.i = data+1;
    b->icon = (data == markers[active_marker].trace) ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    b->bg = LCD_TRACE_1_COLOR+data;
    return;
  }
  markers[active_marker].trace = data;
  menu_move_back(false);
}

static UI_FUNCTION_ADV_CALLBACK(menu_store_trace_acb)
{
  (void)item;
  if(b){
    plot_printf(b->text, sizeof(b->text), S_RARROW"TRACE %d", data+1);
    b->bg = LCD_TRACE_1_COLOR+data;
    if (current_trace == data)
      b->fg = LCD_DARK_GREY;
    return;
  }
  store_trace(current_trace,data);
  menu_move_back(false);
}


static UI_FUNCTION_ADV_CALLBACK(menu_subtract_trace_acb)
{
  (void)item;
  if(b){
    if (data) {
      plot_printf(b->text, sizeof(b->text), "SUBTRACT\nTRACE %d", data);
      b->bg = LCD_TRACE_1_COLOR+data-1;
    }
    else
      plot_printf(b->text, sizeof(b->text), "SUBTRACT\nOFF");
    b->icon = (data == setting.subtract[current_trace]) ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    if (data - 1 == current_trace)
      b->fg = LCD_DARK_GREY;
    return;
  }
  subtract_trace(current_trace,data-1);
  menu_move_back(false);
}

static UI_FUNCTION_ADV_CALLBACK(menu_traces_acb)
{
  (void)item;
  int count = 0;
  if(b){
    if (data == 0) {                // Select trace
      b->param_1.i = current_trace+1;
      b->bg = LCD_TRACE_1_COLOR+current_trace;
    } else if (data == 1) {           // View
      b->icon = IS_TRACE_ENABLE(current_trace) ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    }
    else if (data == 2)               // freeze
      b->icon = setting.stored[current_trace] ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    else if (data == 5) {
      if (setting.subtract[current_trace]) {
        b->icon = BUTTON_ICON_CHECK;
        plot_printf(b->text, sizeof(b->text), "SUBTRACT\nTRACE %d", setting.subtract[current_trace]);
      } else {
        b->icon = BUTTON_ICON_NOCHECK;
        plot_printf(b->text, sizeof(b->text), "SUBTRACT\nOFF");
      }
      // b->icon = setting.subtract[current_trace] ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;  // icon not needed
    } else if (data == 4) {
      if (current_trace == TRACES_MAX-1)
        b->fg = LCD_DARK_GREY;
      else
        b->icon = setting.normalized[current_trace] ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    } else if (data == 3) {
      b->icon = setting.average[current_trace] != AV_OFF ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
      plot_printf(b->text, sizeof(b->text), "CALC\n%s", averageText[setting.average[current_trace]]);
      // b->icon = setting.average[current_trace] ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK; // icon not needed
    }
    return;
  }
  switch(data) {
  case 0:
    menu_push_submenu(menu_trace);
    return;
  case 1:
    for (int i=0;i<TRACES_MAX;i++)
      if (IS_TRACE_ENABLE(i))
        count++;
    if (IS_TRACE_ENABLE(current_trace)) {
      if (count > 1)  {                   // Always 1 trace enabled
        TRACE_DISABLE(1<<current_trace);
      } else {
        drawMessageBox("Trace", "Enable at least one trace", 2000);
        redraw_request|= REDRAW_AREA;
      }
    } else {
      TRACE_ENABLE(1<<current_trace);
    }
    break;
  case 2:                               // Freeze
    setting.stored[current_trace] = !setting.stored[current_trace];
    break;
  case 5:
    menu_push_submenu(menu_subtract_trace);
    return;
    break;
  case 4:
    if (current_trace < TRACES_MAX-1) {
      toggle_normalize(current_trace);
      if (setting.normalized[current_trace]) {
//        kp_help_text = "Ref level";
//        ui_mode_keypad(KM_REFLEVEL);
//        setting.normalize_level = 20.0; // uistat.value;
      } else {
//        set_auto_reflevel(true);
      }
    }
    break;
  case 3:
    menu_push_submenu(menu_average);
    return;
    break;

#ifdef TINYSA4
    case 6:
      save_to_sd(1+(2<<current_trace));      // frequencies + trace
      break;
#endif
  }
//  ui_mode_normal();
//  draw_cal_status();
}
#if 0
static UI_FUNCTION_ADV_CALLBACK(menu_storage_acb)
{
  (void)item;
  if(b){
    if (data == 0 && setting.show_stored)
      b->icon = BUTTON_ICON_CHECK;
    if (setting.subtract[0]){
      if (data == 2 && setting.show_stored)
        b->icon = BUTTON_ICON_CHECK;
      if (data == 3 && !setting.show_stored)
        b->icon = BUTTON_ICON_CHECK;
    }
    return;
  }
  switch(data) {
    case 0:
      store_trace(0,2);
      break;
    case 1:
      set_clear_storage();
      break;
    case 2:
      set_subtract_storage();
      break;
    case 3:
      toggle_normalize();
      if (setting.subtract[0]) {
        kp_help_text = "Ref level";
        ui_mode_keypad(KM_REFLEVEL);
//        setting.normalize_level = uistat.value;
      } else
        set_auto_reflevel(true);
      break;
#ifdef TINYSA4
    case 4:
      save_to_sd(1+2);      // frequencies + actual
      break;
    case 5:
      save_to_sd(1+4);      // frequencies + stored
      break;
#endif
  }
  ui_mode_normal();
//  draw_cal_status();
}
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_waterfall_acb){
  (void)item;
  (void)data;
  if (b){
#ifdef TINYSA4
    if (!(sweep_mode & SWEEP_ENABLE)){
      b->icon = (sweep_mode & SWEEP_ONCE) ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
      plot_printf(b->text, sizeof(b->text), "SINGLE\nSWEEP");
    } else {
      b->icon = setting.waterfall ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
      plot_printf(b->text, sizeof(b->text), "WATER\nFALL");
    }
#else
    b->icon = setting.waterfall ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    plot_printf(b->text, sizeof(b->text), "WATER\nFALL");
#endif
    return;
  }
#ifdef TINYSA4
  if (is_paused()) {
    resume_once(1);
  } else
#endif
  {
    setting.waterfall++; if (setting.waterfall>W_BIG)setting.waterfall = W_OFF;
    if (setting.waterfall != W_OFF)
      setting.level_meter = false;
    set_waterfall();
    ui_mode_normal();
  }
}

#ifdef __LEVEL_METER__
static UI_FUNCTION_ADV_CALLBACK(menu_level_meter_acb){
  (void)item;
  (void)data;
  if (b){
    b->icon = setting.level_meter ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  setting.level_meter = !setting.level_meter;
  if (setting.level_meter)
    setting.waterfall = W_OFF;
  set_level_meter();
  ui_mode_normal();
}
#endif


#ifdef __LIMITS__
uint8_t active_limit = 0;
static UI_FUNCTION_ADV_CALLBACK(menu_limit_select_acb)
{
  (void)item;
  if(b){
    int count = 0;
    for (int i=0;i<LIMITS_MAX;i++) {if (setting.limits[current_trace][i].enabled) count++; }
    if (count == 0) setting.limits[current_trace][0].enabled = true;
    b->icon = (setting.limits[current_trace][data].enabled?BUTTON_ICON_CHECK:BUTTON_ICON_NOCHECK) ;
    plot_printf(b->text, sizeof(b->text), "%.6FHz\n%.2F%s", (float)setting.limits[current_trace][data].frequency, value(setting.limits[current_trace][data].level),unit_string[setting.unit]);
    return;
  }
  active_limit = data;
  setting.limits[current_trace][active_limit].enabled = true;
  dirty = true;
  limits_update();
  menu_push_submenu(menu_limit_modify);
}

#endif

extern const menuitem_t menu_marker_select[];

static UI_FUNCTION_ADV_CALLBACK(menu_marker_modify_acb)
{
  (void)item;
  if (active_marker == MARKER_INVALID) return;
  if(b){
    uistat.text[0] = 0;
    uistat.text[1] = 0;
    switch(data) {
    case 0:
      uistat.text[0] = active_marker+'1';
      break;
    case M_DELTA:
      uistat.text[0] = markers[active_marker].ref+'1';
      /* fall through */
    case M_NOISE:
    case M_TRACKING:
    case M_AVER:
      b->icon = BUTTON_ICON_NOCHECK;
      if (markers[active_marker].mtype & data)
        b->icon = BUTTON_ICON_CHECK;
      break;
    case M_STORED:
      uistat.text[0] = markers[active_marker].trace+'1';
      b->bg = LCD_TRACE_1_COLOR+markers[active_marker].trace;
      break;
    }
    b->param_1.text = uistat.text;
    return;
  }
  if (data == M_DELTA && !(markers[active_marker].mtype & M_DELTA)) {   // Not yet set
    menu_push_submenu(menu_marker_ref_select);
    goto set_delta;
    return;
  } else if (data == M_STORED) {
    current_trace = 0;
    menu_push_submenu(menu_marker_trace);
    return;
  } else if (markers[active_marker].mtype & data)
    markers[active_marker].mtype &= ~data;
  else if (data) {
    set_delta:
    markers[active_marker].mtype |= data;
  } else
    menu_push_submenu(menu_marker_select);
  markmap_all_markers();
//  redraw_marker(active_marker, TRUE);
//  menu_move_back(false);
}

static UI_FUNCTION_ADV_CALLBACK(menu_marker_ref_select_acb)
{
  (void)item;
  if(b){
//    b->icon = markers[data-1].enabled ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    b->icon = (markers[active_marker].ref == data-1 ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP );
    b->param_1.i = data;
    return;
  }
  markers[data-1].enabled = true;
//  interpolate_maximum(data-1);        // possibly not a maximum
  set_marker_index(data-1, markers[data-1].index);
  markers[active_marker].ref = data-1;
  redraw_marker(active_marker);
  menu_move_back(false);
}

extern const menuitem_t menu_marker_modify[];
static UI_FUNCTION_ADV_CALLBACK(menu_marker_select_acb)
{
  (void)item;
  if(b){
    b->icon = markers[data-1].enabled ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    b->param_1.i = data;
    return;
  }
  markers[data-1].enabled = true;
//  interpolate_maximum(data-1);        // possibly not a maximum
  set_marker_index(data-1, markers[data-1].index);
  active_marker_select(data-1);
//  menu_push_submenu(menu_marker_modify);
  redraw_marker(active_marker);
  menu_move_back(false);
}

static UI_FUNCTION_CALLBACK(menu_marker_delete_cb)
{
  (void)item;
  (void)data;
  if (active_marker>=0){
    for (int i = 0; i<MARKER_COUNT; i++ ) {
      if (markers[i].enabled && (markers[i].mtype & M_DELTA) && markers[i].ref == active_marker)
        markers[i].enabled = false;
    }
    markers[active_marker].enabled = false;
    markmap_all_markers();
    menu_move_back(false);
  }
}

#ifdef __LIMITS__
static UI_FUNCTION_CALLBACK(menu_limit_disable_cb)
{
  (void)item;
  (void)data;
  int count = 0;
  for (int i=0;i<LIMITS_MAX;i++) {if (setting.limits[current_trace][i].enabled) count++; }
  if (count == 1 && setting.limits[current_trace][active_limit].enabled) {
    drawMessageBox("Error", "At least one entry",1000);
    return;
  }

  if (active_limit<LIMITS_MAX){
    setting.limits[current_trace][active_limit].enabled = false;
    dirty = true;
    limits_update();
    menu_move_back(false);
  }
}
#endif

#ifdef TINYSA4
static const uint16_t rbwsel_x10[]={0,2,10,30,100,300,1000,3000,6000,8500};
static const char* rbwsel_text[]={"AUTO","200","1k","3k","10k","30k","100k","300k","600k","850k"};
#else
static const uint16_t rbwsel_x10[]={0,30,100,300,1000,3000,6000};
#endif
#ifdef __VBW__
static const uint16_t vbwsel_x100[]={0,100,30,10,3,1};
//static const char* vbwsel_text[]={"AUTO","0.01","0.03", "0.1", "0.3","   "};
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_rbw_acb)
{
  (void)item;
  if (b){
#ifdef TINYSA4
  b->param_1.text = rbwsel_text[data];
#else
    b->param_1.u = rbwsel_x10[data]/10;
#endif
  b->icon = setting.rbw_x10 == rbwsel_x10[data] ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  set_RBW(rbwsel_x10[data]);
  menu_move_back(true);
}

#ifdef __VBW__
static UI_FUNCTION_ADV_CALLBACK(menu_vbw_acb)
{
  (void)item;
  if (b){
  b->param_1.f = vbwsel_x100[data] > 0 ? 1.0f/vbwsel_x100[data] : 0;
  b->icon = setting.vbw_x100 == vbwsel_x100[data] ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  set_VBW(vbwsel_x100[data]);
  menu_move_back(true);
}

#endif

static UI_FUNCTION_ADV_CALLBACK(menu_unit_acb)
{
  (void)item;
  if (b){
    b->icon = data == setting.unit ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  set_unit(data);
  menu_move_back(true);
}

#if 0
enum {
  S_20,S_10,S_5,S_2,S_1,S_P5,S_P2,S_P1,S_P05,S_P02,S_P01
};
static const float menu_scale_per_value[11]={20,10,5,2,1,0.5,0.2,0.1,0.05,0.02,0.01};

static UI_FUNCTION_ADV_CALLBACK(menu_scale_per_acb)
{
  (void)item;
  if(b){
    return;
  }
  set_scale(menu_scale_per_value[data]);
  menu_move_back(true);
}
#endif

const char *mode_text[] = {"PRE","POST","MID"};

static UI_FUNCTION_ADV_CALLBACK(menu_trigger_acb)
{
  (void)item;
  if(b){
    if (data == T_MODE) {
      b->param_1.text = mode_text[setting.trigger_mode - T_PRE];
    } else if (data == T_UP || data == T_DOWN)
      b->icon = setting.trigger_direction == data ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    else
      b->icon = setting.trigger == data ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  if (data == T_MODE) {
    setting.trigger_mode += 1;
    if (setting.trigger_mode > T_MID)
      setting.trigger_mode = T_PRE;
    set_trigger(setting.trigger_mode);
  } else if (data != T_DONE) {
    set_trigger(data);
//  menu_move_back(false);
    ui_mode_normal();
  }
  completed = true;
}

#if 0
static void choose_active_trace(void)
{
  int i;
  if (trace[uistat.current_trace].enabled)
    // do nothing
    return;
  for (i = 0; i < TRACE_COUNT ; i++)
    if (trace[i].enabled) {
      uistat.current_trace = i;
      return;
    }
}
#endif
static void choose_active_marker(void)
{
  int i;
  for (i = 0; i < MARKER_COUNT; i++)
    if (markers[i].enabled == M_ENABLED) {
      active_marker = i;
      return;
    }
  active_marker = MARKER_INVALID;
}
#ifdef TINYSA4
#define __HARMONIC__
#endif

#ifdef __HARMONIC__
static UI_FUNCTION_ADV_CALLBACK(menu_harmonic_acb)
{
  (void)item;
  if(b){
    b->icon = setting.harmonic == data ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  set_harmonic(data);
}
#endif


static UI_FUNCTION_ADV_CALLBACK(menu_settings_agc_acb){
  (void)item;
  (void)data;
  if(b){
    b->icon = AUTO_ICON(setting.agc);
    return;
  }
  toggle_AGC();
}

static UI_FUNCTION_ADV_CALLBACK(menu_settings_lna_acb){
  (void)item;
  (void)data;
  if(b){
    if (S_IS_AUTO(setting.lna))
      b->icon = BUTTON_ICON_CHECK_AUTO;
    else
      b->icon = setting.lna ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_LNA();
}

static UI_FUNCTION_ADV_CALLBACK(menu_settings_bpf_acb){
  (void)item;
  (void)data;
  if(b){
    b->icon = setting.tracking ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_tracking();
}

#ifdef __LCD_BRIGHTNESS__
static UI_FUNCTION_CALLBACK(menu_brightness_cb)
{
  (void)item;
  (void)data;
  int16_t value = config._brightness;
  ili9341_set_foreground(LCD_MENU_TEXT_COLOR);
  ili9341_set_background(LCD_MENU_COLOR);
  ili9341_fill(LCD_WIDTH/2-80, LCD_HEIGHT/2-20, 160, 40);
  ili9341_drawstring_7x13("BRIGHTNESS", LCD_WIDTH/2-35, LCD_HEIGHT/2-13);
  ili9341_drawstring_7x13(S_LARROW" USE LEVELER BUTTON "S_RARROW, LCD_WIDTH/2-72, LCD_HEIGHT/2+2);
  while (TRUE) {
    int status = btn_check();
    if (status & (EVT_UP|EVT_DOWN)) {
      do {
        if (status & EVT_UP  ) value+=5;
        if (status & EVT_DOWN) value-=5;
        if (value <   0) value =   0;
        if (value > 100) value = 100;
        lcd_setBrightness(value);
        status = btn_wait_release();
      } while (status != 0);
    }
    if (status == EVT_BUTTON_SINGLE_CLICK)
      break;
  }
  config._brightness = (uint8_t)value;
  lcd_setBrightness(value);
  redraw_request|= REDRAW_AREA;
  ui_mode_normal();
}
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_settings_pulse_acb){
  (void)item;
  (void)data;
  if(b){
    b->icon = setting.pulse ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_pulse();
}

#ifdef __DRAW_LINE__
static UI_FUNCTION_ADV_CALLBACK(menu_settings_draw_line_acb){
  (void)item;
  (void)data;
  if(b){
    b->icon = setting.draw_line ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_draw_line();
  if (setting.draw_line) {
    kp_help_text = "Level";
    ui_mode_keypad(KM_TRIGGER);
    set_trigger(T_AUTO);
  }
}
#endif
#ifdef __HAM_BAND__
static UI_FUNCTION_ADV_CALLBACK(menu_settings_ham_bands){
  (void)item;
  (void)data;
  if(b){
    b->icon = config.hambands ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_hambands();
}
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_settings_below_if_acb){
  (void)item;
  (void)data;
  if(b){
    if (S_IS_AUTO(setting.below_IF))
      b->icon = BUTTON_ICON_CHECK_AUTO;
    else
      b->icon = setting.below_IF ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_below_IF();
}

#if 0
static UI_FUNCTION_ADV_CALLBACK(menu_settings_ultra_acb){
  (void)item;
  (void)data;
  if(b){
    if (S_IS_AUTO(setting.ultra))
      b->icon = BUTTON_ICON_CHECK_AUTO;
    else
      b->icon = setting.ultra ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_ultra();
}
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_lo_output_acb){
  (void)item;
  (void)data;
  if(b){
    b->icon = setting.tracking_output ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_tracking_output();
}

static UI_FUNCTION_ADV_CALLBACK(menu_pause_acb)
{
  (void) data;
  (void) item;
  if (b){
    b->icon = !(sweep_mode & SWEEP_ENABLE) ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  toggle_sweep();
//  menu_move_back(true);
//  draw_cal_status();
}

static UI_FUNCTION_ADV_CALLBACK(menu_flip_acb)
{
  (void) data;
  (void) item;
  if (b){
    return;
  }
  config.flip = ! config.flip;
  ili9341_flip(config.flip);
  config_save();
  redraw_request|= REDRAW_AREA | REDRAW_FREQUENCY | REDRAW_CAL_STATUS;
}

static UI_FUNCTION_ADV_CALLBACK(menu_shift_acb)
{
  (void) data;
  (void) item;
  if (b){
    b->icon = setting.frequency_offset != FREQUENCY_SHIFT ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  if (setting.frequency_offset != FREQUENCY_SHIFT) {
    setting.frequency_offset = FREQUENCY_SHIFT;
  } else {
    if (FREQ_IS_STARTSTOP()) {
      freq_t old_start = get_sweep_frequency(ST_START);
      freq_t old_stop = get_sweep_frequency(ST_STOP);
      kp_help_text = "Actual start frequency";
      ui_mode_keypad(KM_START);
      setting.frequency_offset = uistat.value - old_start + FREQUENCY_SHIFT;
      set_sweep_frequency(ST_START, old_start);
      set_sweep_frequency(ST_STOP, old_stop);
    } else {
      freq_t old_center = get_sweep_frequency(ST_CENTER);
      freq_t old_span = get_sweep_frequency(ST_SPAN);
      kp_help_text = "Actual center frequency";
      ui_mode_keypad(KM_CENTER);
      setting.frequency_offset = uistat.value - old_center + FREQUENCY_SHIFT;
      set_sweep_frequency(ST_CENTER, old_center);
      set_sweep_frequency(ST_SPAN, old_span);
    }
  }
  ui_mode_normal();
//  menu_move_back(true);
//  draw_cal_status();
}

#ifdef __REMOTE_DESKTOP__
#if 0   // Not used in UI
static UI_FUNCTION_ADV_CALLBACK(menu_send_display_acb)
{
  (void) data;
  (void) item;
  if (b){
    b->icon = auto_capture ? BUTTON_ICON_CHECK : BUTTON_ICON_NOCHECK;
    return;
  }
  auto_capture = ! auto_capture;
  // Update all screen to CPU
  if (auto_capture)
    redraw_request|=REDRAW_AREA|REDRAW_BATTERY|REDRAW_FREQUENCY|REDRAW_CAL_STATUS;
}
#endif
#endif

static UI_FUNCTION_ADV_CALLBACK(menu_outputmode_acb)
{
  (void) data;
  (void) item;
  if(b){
    b->param_1.text = setting.mute ? "OFF" : "ON";
    return;
  }
  toggle_mute();
}

static UI_FUNCTION_ADV_CALLBACK(menu_enter_marker_acb)
{
  (void) data;
  (void) item;
  if(b){
    b->param_1.text = FREQ_IS_CW() ? "TIME" : "FREQUENCY";
    return;
  }
  if (FREQ_IS_CW())
    ui_mode_keypad(KM_MARKER_TIME);
  else
    ui_mode_keypad(KM_MARKER);
}

#ifdef TINYSA4
static const uint16_t points_setting[] = {51, 101, 201, 256, 290, 450};
#else
static const uint16_t points_setting[] = {51, 101, 145, 290};
#endif
static UI_FUNCTION_ADV_CALLBACK(menu_points_acb){
  (void)item;
  if(b){
    b->icon = points_setting[data] == sweep_points ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    b->param_1.i = points_setting[data];
    return;
  }
  set_sweep_points(points_setting[data]);
}

#ifdef __USE_SERIAL_CONSOLE__
static UI_FUNCTION_ADV_CALLBACK(menu_serial_speed_acb)
{
  static const uint32_t usart_speed[] = {19200, 38400, 57600, 115200, 230400, 460800, 921600, 1843200, 2000000, 3000000};
  (void)item;
  uint32_t speed = usart_speed[data];
  if (b){
    b->icon = config._serial_speed == speed ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    b->param_1.u = speed;
    return;
  }
  config._serial_speed = speed;
  shell_update_speed();
}

static UI_FUNCTION_ADV_CALLBACK(menu_connection_acb)
{
  (void)item;
  if (b){
    b->icon = (config._mode&_MODE_CONNECTION_MASK) == data ? BUTTON_ICON_GROUP_CHECKED : BUTTON_ICON_GROUP;
    return;
  }
  config._mode&=~_MODE_CONNECTION_MASK;
  config._mode|=data;
  config_save();
  shell_reset_console();
}
#endif
// ===[MENU DEFINITION]=========================================================
// Back button submenu list

static const menuitem_t menu_back[] = {
  { MT_CANCEL,   0, S_LARROW" BACK", NULL },
  { MT_NONE,     0, NULL, NULL } // sentinel
};

static const menuitem_t menu_store_preset[] =
{
  { MT_ADV_CALLBACK, 0,  "STORE AS\nSTARTUP",menu_store_preset_acb},
  { MT_ADV_CALLBACK |MT_REPEATS,  DATA_STARTS_REPEATS(1,4),  "STORE %d",         menu_store_preset_acb},
  { MT_ADV_CALLBACK, 100,"FACTORY\nDEFAULTS",menu_store_preset_acb},
  { MT_NONE,     0,     NULL,menu_back} // next-> menu_back
};

static const menuitem_t menu_load_preset[] =
{
  { MT_ADV_CALLBACK,            0,                          "LOAD\nSTARTUP", menu_load_preset_acb},
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(1,4),   MT_CUSTOM_LABEL, menu_load_preset_acb},
  { MT_SUBMENU,                 0,                          "STORE"  ,       menu_store_preset},
  { MT_NONE,     0,     NULL, menu_back} // next-> menu_back
};
#ifdef TINYSA4
static const menuitem_t menu_mixer_drive[] = {
  { MT_ADV_CALLBACK, 4, "Auto",     menu_mixer_drive_acb},
  { MT_ADV_CALLBACK, 3, "%+ddBm",   menu_mixer_drive_acb},
  { MT_ADV_CALLBACK, 2, "%+ddBm",   menu_mixer_drive_acb},
  { MT_ADV_CALLBACK, 1, "%+ddBm",   menu_mixer_drive_acb},
  { MT_ADV_CALLBACK, 0, "%+ddBm",   menu_mixer_drive_acb},
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};
#else
static const menuitem_t menu_lo_drive[] = {
  { MT_ADV_CALLBACK, 15, "%+ddBm",   menu_lo_drive_acb},
  { MT_ADV_CALLBACK, 14, "%+ddBm",   menu_lo_drive_acb},
  { MT_ADV_CALLBACK, 13, "%+ddBm",   menu_lo_drive_acb},
  { MT_ADV_CALLBACK, 12, "%+ddBm",   menu_lo_drive_acb},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};
#endif

static const menuitem_t  menu_modulation[] = {
  { MT_FORM | MT_TITLE,    0,  "MODULATION",NULL},
  { MT_FORM | MT_ADV_CALLBACK, MO_NONE,             MT_CUSTOM_LABEL,    menu_modulation_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_LOW, MO_AM,      "AM",    menu_modulation_acb},
#ifdef TINYSA4
  { MT_FORM | MT_ADV_CALLBACK, MO_WFM,              "FM",               menu_modulation_acb},
  { MT_FORM | MT_KEYPAD,   KM_MODULATION,           "FREQ: %s",         "50Hz..5kHz"},
  { MT_FORM | MT_KEYPAD,   KM_DEPTH,               "AM DEPTH: %s%%",         "0..100"},
  { MT_FORM | MT_KEYPAD,   KM_DEVIATION,            "FM DEVIATION: %s",         "1kHz..300kHz"},
//  { MT_FORM | MT_ADV_CALLBACK, MO_NFM2,              MT_CUSTOM_LABEL,    menu_modulation_acb},
//  { MT_FORM | MT_ADV_CALLBACK, MO_NFM3,              MT_CUSTOM_LABEL,    menu_modulation_acb},
#else
  { MT_FORM | MT_ADV_CALLBACK, MO_NFM,              MT_CUSTOM_LABEL,    menu_modulation_acb},
  { MT_FORM | MT_ADV_CALLBACK, MO_WFM,              MT_CUSTOM_LABEL,    menu_modulation_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_LOW, MO_EXTERNAL,MT_CUSTOM_LABEL,    menu_modulation_acb},
  { MT_FORM | MT_KEYPAD,   KM_MODULATION,           "FREQ: %s",         "50Hz..5kHz"},
#endif
  { MT_FORM | MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t  menu_sweep[] = {
  { MT_FORM | MT_KEYPAD,   KM_SPAN,             "SPAN: %s",         VARIANT("0..350MHz", range_text)},
  { MT_FORM | MT_KEYPAD | MT_LOW, KM_LEVELSWEEP,"LEVEL CHANGE: %s", VARIANT("-70..70","-90..90")},
  { MT_FORM | MT_KEYPAD,   KM_SWEEP_TIME,       "SWEEP TIME: %s",   "0..600 seconds"},
  { MT_FORM | MT_SUBMENU,  0,                   "SWEEP POINTS",   menu_sweep_points_form},
  { MT_FORM | MT_NONE, 0, NULL, menu_back} // next-> menu_back
};
#ifdef TINYSA4
static const menuitem_t  menu_lowoutput_settings[] = {
  { MT_FORM | MT_ADV_CALLBACK,  0,              "Cleanest signal, max 4.4GHz",       menu_lowoutput_settings_acb},
  { MT_FORM | MT_ADV_CALLBACK,  1,              "Highest accuracy, max 5.4GHz",    menu_lowoutput_settings_acb},
  { MT_FORM | MT_SUBMENU,  255, S_RARROW"Config", menu_config},
  { MT_FORM | MT_NONE, 0, NULL, menu_back} // next-> menu_back
};
#endif

char low_level_help_text[12] = "-76..-6";
char center_text[18] = "FREQ: %s";

const menuitem_t  menu_lowoutputmode[] = {
  { MT_FORM | MT_ADV_CALLBACK, 0,               "LOW OUTPUT            %s", menu_outputmode_acb},
//  { MT_FORM | MT_ADV_CALLBACK,  0,              "MOD: %s",   menu_smodulation_acb},
  { MT_FORM | MT_KEYPAD,   KM_CENTER,           center_text,         VARIANT("10kHz..350MHz","10kHz..850MHz")},
  { MT_FORM | MT_KEYPAD,   KM_LOWOUTLEVEL,      "LEVEL: %s",        low_level_help_text},
  { MT_FORM | MT_ADV_CALLBACK,  0,              MT_CUSTOM_LABEL,   menu_smodulation_acb},
  { MT_FORM | MT_ADV_CALLBACK,  0,              MT_CUSTOM_LABEL,   menu_sweep_acb},
#ifdef __SWEEP_RESTART__
  { MT_FORM | MT_ADV_CALLBACK,  0,              MT_CUSTOM_LABEL,   menu_restart_acb},
#endif
  { MT_FORM | MT_KEYPAD,  KM_EXT_GAIN,            "EXTERNAL GAIN: %s",   "-100..+100"},
#ifdef TINYSA4
  { MT_FORM | MT_ADV_CALLBACK,  255, "OUTPUT: %s",     menu_lowoutput_settings_acb},
#endif
  { MT_FORM | MT_CANCEL,   0,                   "MODE",             NULL },
  { MT_FORM | MT_NONE, 0, NULL, NULL } // sentinel
};

const menuitem_t  menu_highoutputmode[] = {
  { MT_FORM | MT_ADV_CALLBACK,  0,      "HIGH OUTPUT           %s", menu_outputmode_acb},
  { MT_FORM | MT_KEYPAD,    KM_CENTER,  center_text,         VARIANT("240MHz..959MHz",range_text)},
  { MT_FORM | MT_KEYPAD,   KM_HIGHOUTLEVEL,      "LEVEL: %s",        low_level_help_text /* "-76..-6" */},
  { MT_FORM | MT_ADV_CALLBACK,   0,     MT_CUSTOM_LABEL,   menu_smodulation_acb},
  { MT_FORM | MT_ADV_CALLBACK,  0,      MT_CUSTOM_LABEL,   menu_sweep_acb},
#ifdef __SWEEP_RESTART__
  { MT_FORM | MT_ADV_CALLBACK,  0,      MT_CUSTOM_LABEL,   menu_restart_acb},
#endif
  { MT_FORM | MT_KEYPAD,  KM_EXT_GAIN,            "EXTERNAL GAIN: %s",          "-100..+100"},
#ifdef TINYSA4
  { MT_FORM | MT_SUBMENU,  255, S_RARROW" Settings", menu_settings3},
#endif
  { MT_FORM | MT_CANCEL,    0,          "MODE",             NULL },
  { MT_FORM | MT_NONE, 0, NULL, NULL } // sentinel
};

static const menuitem_t  menu_average[] = {
  { MT_ADV_CALLBACK, AV_OFF,        "OFF",          menu_average_acb},
  { MT_ADV_CALLBACK, AV_MIN,        "MIN\nHOLD",    menu_average_acb},
  { MT_ADV_CALLBACK, AV_MAX_HOLD,   "MAX\nHOLD",    menu_average_acb},
  { MT_ADV_CALLBACK, AV_MAX_DECAY,  "MAX\nDECAY",   menu_average_acb},
  { MT_ADV_CALLBACK, AV_4,          "AVER 4",       menu_average_acb},
  { MT_ADV_CALLBACK, AV_16,         "AVER 16",      menu_average_acb},
#ifdef TINYSA4
  { MT_ADV_CALLBACK, AV_100,        "AVER",     menu_average_acb},
#endif
#ifdef __QUASI_PEAK__
  { MT_ADV_CALLBACK, AV_QUASI,      "QUASI\nPEAK",  menu_average_acb},
#endif
#ifdef __LIMITS__
  { MT_ADV_CALLBACK, AV_TABLE,      "TABLE"S_RARROW"\nTRACE",  menu_average_acb},
#endif
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_rbw[] = {
  { MT_ADV_CALLBACK, 0, "  AUTO",   menu_rbw_acb},
#ifdef TINYSA4
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(1,9), "%sHz",   menu_rbw_acb},
#else
  { MT_ADV_CALLBACK | MT_REPEATS, DATA_STARTS_REPEATS(1,6), "%4dkHz",   menu_rbw_acb},
#endif
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back

};

#ifdef __VBW__
static const menuitem_t menu_vbw[] = {
  { MT_ADV_CALLBACK, 0, "     AUTO",   menu_vbw_acb},
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(1,5), "%b.2f RBW",   menu_vbw_acb},
  { MT_NONE,      0, NULL, menu_back} // next-> menu_back
};
#endif

static const menuitem_t menu_reffer[] = {
  { MT_FORM | MT_ADV_CALLBACK|MT_REPEATS,  DATA_STARTS_REPEATS(0,8), "%s", menu_reffer_acb},
  { MT_FORM | MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_atten[] = {
  { MT_ADV_CALLBACK,          0,           "%s",    menu_atten_acb},
  { MT_KEYPAD | MT_LOW,   KM_ATTENUATION,  "MANUAL\n\b%s",  "0 - 30dB"},
//  { MT_ADV_CALLBACK | MT_HIGH,0,           "0dB",     menu_atten_high_acb},
  { MT_ADV_CALLBACK | MT_HIGH,30,          "22.5 - 40dB",    menu_atten_high_acb},
  { MT_FORM | MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_reflevel[] = {
  { MT_ADV_CALLBACK,0,        "AUTO",    menu_reflevel_acb},
  { MT_KEYPAD,  KM_REFLEVEL,  "MANUAL\n\b%s",  NULL},
  { MT_CANCEL, 0,      S_LARROW" BACK", NULL },
  { MT_NONE,   0, NULL, NULL } // sentinel
};

const menuitem_t menu_marker_search[] = {
  { MT_CALLBACK, 4, "PEAK\n SEARCH",          menu_marker_search_cb },
  { MT_CALLBACK, 0, "MIN\n" S_LARROW" LEFT",  menu_marker_search_cb },
  { MT_CALLBACK, 1, "MIN\n" S_RARROW" RIGHT", menu_marker_search_cb },
  { MT_CALLBACK, 2, "MAX\n" S_LARROW" LEFT",  menu_marker_search_cb },
  { MT_CALLBACK, 3, "MAX\n" S_RARROW" RIGHT", menu_marker_search_cb },
  { MT_ADV_CALLBACK, 0, "ENTER\n%s",          menu_enter_marker_acb},
  { MT_ADV_CALLBACK, M_TRACKING,    "TRACKING",menu_marker_modify_acb },
#ifdef TINYSA4
  { MT_KEYPAD,   KM_NOISE,      "PEAK\n\b%s",   "2..20 dB"},
#endif
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

const menuitem_t menu_marker_modify[] = {
  { MT_ADV_CALLBACK, 0,            "MARKER %s",      menu_marker_modify_acb},
  { MT_ADV_CALLBACK, M_DELTA,       "DELTA %s",    menu_marker_modify_acb},
  { MT_ADV_CALLBACK, M_NOISE,       "NOISE",    menu_marker_modify_acb},
  { MT_ADV_CALLBACK, M_TRACKING,    "TRACKING", menu_marker_modify_acb},
  { MT_ADV_CALLBACK, M_STORED,      "TRACE %s",   menu_marker_modify_acb},
  { MT_ADV_CALLBACK, M_AVER,      "TRACE\nAVERAGE",   menu_marker_modify_acb},
  { MT_SUBMENU,  0,                 "SEARCH",   menu_marker_search},
  { MT_CALLBACK, M_DELETE,          "DELETE",   menu_marker_delete_cb},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

#ifdef __LIMITS__
static const menuitem_t menu_limit_modify[] =
{
  { MT_KEYPAD,   KM_LIMIT_FREQ,   "FREQUENCY\n\b%s",          "Frequency"},
  { MT_KEYPAD,   KM_LIMIT_LEVEL,  "LEVEL\n\b%s",              "Level"},
  { MT_CALLBACK, 0,               "DISABLE",            menu_limit_disable_cb},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_limit_select[] = {
  { MT_ADV_CALLBACK | MT_REPEATS,   DATA_STARTS_REPEATS(0,LIMITS_MAX), MT_CUSTOM_LABEL, menu_limit_select_acb },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};
#endif

#if 0
const menuitem_t menu_marker_sel[] = {
  { MT_CALLBACK, 1, "MARKER %d", menu_marker_sel_cb },
  { MT_CALLBACK, 2, "MARKER %d", menu_marker_sel_cb },
  { MT_CALLBACK, 3, "MARKER %d", menu_marker_sel_cb },
  { MT_CALLBACK, 4, "MARKER %d", menu_marker_sel_cb },
//  { MT_CALLBACK, 0, "ALL OFF", menu_marker_sel_cb },
  { MT_CALLBACK, 0, "DELTA", menu_marker_sel_cb },
  { MT_CALLBACK, 0, "NOISE", menu_marker_sel_cb },
  { MT_CALLBACK, 0, "TRACKING", menu_marker_sel_cb },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};
#endif

const menuitem_t menu_marker_select[] = {
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(1,MARKER_COUNT), "MARKER %d", menu_marker_select_acb },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_marker_ref_select[] = {
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(1,MARKER_COUNT), S_RARROW"MARKER %d", menu_marker_ref_select_acb },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

const menuitem_t menu_marker_ops[] = {
  { MT_CALLBACK, ST_START,  S_RARROW" START",    menu_marker_op_cb },
  { MT_CALLBACK, ST_STOP,   S_RARROW" STOP",     menu_marker_op_cb },
  { MT_CALLBACK, ST_CENTER, S_RARROW" CENTER",   menu_marker_op_cb },
  { MT_CALLBACK, ST_SPAN,   S_RARROW" SPAN",     menu_marker_op_cb },
  { MT_CALLBACK, 4,         S_RARROW" REF LEVEL",menu_marker_op_cb },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_marker[] = {
//  { MT_SUBMENU,  0, "SELECT\nMARKER",     menu_marker_sel},
  { MT_SUBMENU,  0, "MODIFY\nMARKERS",    menu_marker_modify},
  { MT_SUBMENU,  0, "MARKER\nOPS", menu_marker_ops},
  { MT_SUBMENU,  0, "SEARCH\nMARKER",     menu_marker_search},
  { MT_CALLBACK, 0, "RESET\nMARKERS",     menu_markers_reset_cb},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

#ifndef TINYSA4
static const menuitem_t menu_dfu[] = {
  { MT_FORM | MT_CALLBACK, 0, "ENTER DFU",      menu_dfu_cb},
  { MT_FORM | MT_NONE,     0, NULL, menu_back} // next-> menu_back
};
#endif

#ifdef __HARMONIC__
static const menuitem_t menu_harmonic[] =
{
  { MT_ADV_CALLBACK, 0,   "OFF",              menu_harmonic_acb},
  { MT_ADV_CALLBACK, 2,     "2",              menu_harmonic_acb},
  { MT_ADV_CALLBACK, 3,     "3",              menu_harmonic_acb},
  { MT_ADV_CALLBACK, 4,     "4",              menu_harmonic_acb},
  { MT_ADV_CALLBACK, 5,     "5",              menu_harmonic_acb},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};
#endif


static const menuitem_t menu_scanning_speed[] =
{
// { MT_ADV_CALLBACK,     SD_NORMAL,  "NORMAL",           menu_scanning_speed_acb},    // order must match definition of enum
// { MT_ADV_CALLBACK,     SD_PRECISE, PRECISE",           menu_scanning_speed_acb},
// { MT_ADV_CALLBACK | MT_LOW,SD_FAST,  "FAST",           menu_scanning_speed_acb},
// { MT_KEYPAD   | MT_LOW,KM_FAST_SPEEDUP,    "FAST\nSPEEDUP",   "2..20"},
  { MT_KEYPAD, KM_SAMPLETIME,        "SDELAY\n\b%s",   "250..10000, 0=auto"},              // This must be item 4 to match highlighting
  { MT_KEYPAD, KM_OFFSET_DELAY,      "ODELAY\n\b%s",   "250..10000, 0=auto"},              // This must be item 5 to match highlighting
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_sweep_points[] = {
#ifdef TINYSA4
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(0,6), "%3d point", menu_points_acb },
#else
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(0,4), "%3d point", menu_points_acb },
#endif
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_sweep_points_form[] = {
#ifdef TINYSA4
  { MT_FORM|MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(0,6), "%3d point", menu_points_acb },
#else
  { MT_FORM|MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(0,4), "%3d point", menu_points_acb },
#endif
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_sweep_speed[] =
{
 { MT_ADV_CALLBACK,     SD_NORMAL,     "NORMAL",          menu_scanning_speed_acb},    // order must match definition of enum
 { MT_ADV_CALLBACK,     SD_PRECISE,    "PRECISE",         menu_scanning_speed_acb},
#ifdef TINYSA4
 { MT_ADV_CALLBACK,     SD_FAST,        "FAST",            menu_scanning_speed_acb},
#else
 { MT_ADV_CALLBACK | MT_LOW,SD_FAST,   "FAST",            menu_scanning_speed_acb},
#endif
 { MT_ADV_CALLBACK,     SD_NOISE_SOURCE,      "NOISE\nSOURCE",   menu_scanning_speed_acb},
#ifdef TINYSA4
 { MT_KEYPAD,           KM_FAST_SPEEDUP,"SPEEDUP\n\b%s",  "2..20, 0=disable"},
#else
 { MT_KEYPAD   | MT_LOW,KM_FAST_SPEEDUP,"SPEEDUP\n\b%s",  "2..20, 0=disable"},
#endif
 { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

#ifdef TINYSA4
static const menuitem_t menu_curve3[] = {
  { MT_FORM | MT_ADV_CALLBACK | MT_REPEATS, DATA_STARTS_REPEATS(14,6), MT_CUSTOM_LABEL, menu_curve_acb },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_curve2[] = {
  { MT_FORM | MT_ADV_CALLBACK | MT_REPEATS, DATA_STARTS_REPEATS(7,7), MT_CUSTOM_LABEL, menu_curve_acb },
  { MT_FORM | MT_SUBMENU,      0,  S_RARROW" MORE",     menu_curve3},
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_curve[] = {
  { MT_FORM | MT_ADV_CALLBACK | MT_REPEATS, DATA_STARTS_REPEATS(0,7), MT_CUSTOM_LABEL, menu_curve_acb },
  { MT_FORM | MT_SUBMENU,      0,  S_RARROW" MORE",     menu_curve2},
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_curve_confirm[] = {
  { MT_CALLBACK, 1,               "OK",       menu_curve_confirm_cb },
  { MT_CALLBACK, 0,               "CANCEL",   menu_curve_confirm_cb },
  { MT_NONE, 0, NULL, NULL } // sentinel
};
#if 0
static const menuitem_t menu_noise_figure_confirm[] = {
  { MT_CALLBACK, 1,               "STORE\nTINYSA NF",       menu_noise_figure_confirm_cb },
  { MT_CALLBACK, 0,               "CANCEL",   menu_noise_figure_confirm_cb },
  { MT_NONE, 0, NULL, NULL } // sentinel
};
#endif
#endif

#ifdef TINYSA4
static const menuitem_t menu_actual_power2[] =
{
 { MT_ADV_CALLBACK,     0,              "30MHz\nLEVEL", menu_output_level_acb},
 { MT_ADV_CALLBACK,     0,              "1GHz\nLEVEL", menu_output_level2_acb},
 { MT_ADV_CALLBACK,     0,              "1.2GHz\nLEVEL", menu_output_level3_acb},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};
#endif

static const menuitem_t menu_actual_power[] =
{
 { MT_KEYPAD,           KM_ACTUALPOWER, "INPUT\nLEVEL",  "Enter actual level under marker"},
#ifdef TINYSA4
 { MT_SUBMENU,      0,                  "OUTPUT\nLEVEL", menu_actual_power2},
 { MT_CALLBACK,     0,                  "INPUT\nCURVE",  menu_input_curve_prepare_cb},
 { MT_CALLBACK,     0,                  "LNA\nCURVE",    menu_lna_curve_prepare_cb},
 { MT_CALLBACK,     0,                  "ULTRA\nCURVE",  menu_ultra_curve_prepare_cb},
 { MT_CALLBACK,     0,                  "LNA_U\nCURVE",    menu_lna_u_curve_prepare_cb},
 { MT_CALLBACK,     0,                  "OUTPUT\nCURVE", menu_output_curve_prepare_cb},
#else
 { MT_ADV_CALLBACK,     0,              "OUTPUT\nLEVEL", menu_output_level_acb},
#endif
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

#ifdef TINYSA4
static const menuitem_t menu_settings4[];
#endif

#ifdef TINYSA4
static const menuitem_t menu_settings4[] =
{
 { MT_ADV_CALLBACK,     0,     "DEBUG\nFREQ",          menu_debug_freq_acb},
 { MT_ADV_CALLBACK,     0,     "DEBUG\nAVOID",          menu_debug_avoid_acb},
 { MT_ADV_CALLBACK,     0,     "DEBUG\nSPUR",        menu_debug_spur_acb},
 { MT_ADV_CALLBACK,     0,     "HIDE\n21MHz",        menu_hide_21MHz_acb},
#if 0                                                                           // only used during development
  { MT_KEYPAD,   KM_COR_AM,     "COR\nAM", "Enter AM modulation correction"},
  { MT_KEYPAD,   KM_COR_WFM,     "COR\nWFM", "Enter WFM modulation correction"},
  { MT_KEYPAD,   KM_COR_NFM,     "COR\nNFM", "Enter NFM modulation correction"},
#endif
//  { MT_CALLBACK,        0 ,     "CLEAR\nCONFIG",    menu_clearconfig_cb},
//  { MT_SUBMENU,  0,             S_RARROW" MORE",     menu_settings3},
  { MT_KEYPAD,   KM_DIRECT_START,     "DSTART\n\b%s", ""},
  { MT_KEYPAD,   KM_DIRECT_STOP,     "DSTOP\n\b%s", ""},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};
#endif

static const menuitem_t menu_settings3[] =
{
#ifdef TINYSA4
//  { MT_KEYPAD,   KM_GRIDLINES,  "MINIMUM\nGRIDLINES", "Enter minimum horizontal grid divisions"},
#ifndef __NEW_SWITCHES__
  { MT_ADV_CALLBACK,     0,     "ADF OUT",          menu_adf_out_acb},
#endif
  { MT_KEYPAD,   KM_ULTRA_START,"ULTRASTART\n\b%s",   "10G=auto"},
//  { MT_KEYPAD | MT_LOW, KM_IF2,  "IF2 FREQ",           "Set to zero for no IF2"},
  { MT_KEYPAD,  KM_R,  "R\n\b%s",           "Set R"},
  { MT_KEYPAD,  KM_MOD,  "MODULO\n\b%s",           "Set MODULO"},
  { MT_KEYPAD,  KM_CP,  "CP",             "Set CP"},
#ifdef __WAIT_CTS_WHILE_SLEEPING__
  { MT_ADV_CALLBACK,     0,     "SLEEP\nWAIT",    menu_sleep_acb},
#endif
#ifdef __HARMONIC__
  { MT_SUBMENU          ,0,               "HARMONIC",         menu_harmonic},
#endif
//  { MT_ADV_CALLBACK | MT_LOW, 0,    "ULTRA\nMODE",      menu_settings_ultra_acb},
#ifdef __HAM_BAND__
  { MT_ADV_CALLBACK, 0,         "HAM\nBANDS",         menu_settings_ham_bands},
#endif
  { MT_SUBMENU,  0,             S_RARROW" MORE",     menu_settings4},
#else
#ifdef __ULTRA__
  { MT_ADV_CALLBACK,     0,     "ENABLE\nULTRA",    menu_ultra_acb},
  { MT_KEYPAD,   KM_ULTRA_START,        "ULTRA\nSTART",   "Enter ULTRA mode start freq"},
  { MT_ADV_CALLBACK,     0,     "DEBUG\nSPUR",        menu_debug_spur_acb},
#endif
#ifdef __HARMONIC__
  { MT_SUBMENU | MT_HIGH,0,               "HARMONIC",         menu_harmonic},
//  { MT_ADV_CALLBACK,0,          "SPUR\nREMOVAL",          menu_harmonic_spur_acb},
#endif
#ifdef __HAM_BAND__
  { MT_ADV_CALLBACK, 0,         "HAM\nBANDS",         menu_settings_ham_bands},
#endif
#endif  // TINYSA4
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_settings2[] =
{
  { MT_ADV_CALLBACK, 0,             "AGC",           menu_settings_agc_acb},
  { MT_ADV_CALLBACK, 0,             "LNA",           menu_settings_lna_acb},
  { MT_ADV_CALLBACK | MT_LOW, 0,    "BPF",           menu_settings_bpf_acb},
  { MT_ADV_CALLBACK | MT_LOW, 0,    "BELOW IF",      menu_settings_below_if_acb},
  { MT_KEYPAD | MT_LOW, KM_IF,  "IF FREQ\n\b%s",           "0=auto IF"},
#ifdef TINYSA4
  #ifdef __QUASI_PEAK__
  { MT_KEYPAD,   KM_DECAY,      "DECAY\n\b%s",   "0..1000000ms or sweeps"},
  { MT_KEYPAD,   KM_ATTACK,      "ATTACK\n\b%s",   "0..100000ms"},
#endif
#endif
  { MT_SUBMENU,0,               "SCAN\nSPEED",        menu_scanning_speed},
#ifdef TINYSA4
  { MT_SUBMENU | MT_LOW,0,      "MIXER\nDRIVE",      menu_mixer_drive},
  { MT_SUBMENU,  0,             S_RARROW" MORE",     menu_settings3},
#else
  { MT_SUBMENU | MT_LOW,0,      "MIXER\nDRIVE",      menu_lo_drive},
  { MT_KEYPAD,   KM_10MHZ,      "CORRECT\nFREQUENCY", "Enter actual l0MHz frequency"},
#endif
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};


static const menuitem_t menu_settings[] =
{
#ifdef TINYSA4
  { MT_ADV_CALLBACK,0,              "PROGRESS\nBAR",        menu_progress_bar_acb},
  { MT_ADV_CALLBACK,     0,         "DIRECT\nMODE",         menu_direct_acb},
  { MT_ADV_CALLBACK,     0,         "LINEAR\nAVERAGING",    menu_linear_averaging_acb},
  { MT_KEYPAD,      KM_FREQ_CORR,   "FREQ CORR\n\b%s",      "Enter ppb correction"},
//  { MT_SUBMENU,     0,              "CALIBRATE\nHARMONIC",  menu_calibrate_harmonic},
#endif
#ifdef __NOISE_FIGURE__
  { MT_KEYPAD,      KM_NF,          "NF\n\b%s",             "Enter tinySA noise figure"},
#endif
#ifdef __SD_CARD_LOAD__
  { MT_CALLBACK,    0 ,             "LOAD\nCONFIG.INI",     menu_load_config_cb},
//  { MT_CALLBACK,        1 ,       "LOAD\nSETTING.INI",    menu_load_config_cb},
#endif
#ifdef TINYSA4
  { MT_ADV_CALLBACK,     0,              "INTERNALS",            menu_internals_acb},
#endif
  { MT_NONE,        0, NULL, menu_back} // next-> menu_back
};

#ifdef __NOISE_FIGURE__
static const menuitem_t menu_measure_noise_figure[] =
{
 { MT_ADV_CALLBACK,            M_NF_TINYSA,         "MEASURE\nTINYSA NF",menu_measure_acb},
 { MT_ADV_CALLBACK,            M_NF_STORE,          "STORE\nTINYSA NF",menu_measure_acb},
 { MT_ADV_CALLBACK,            M_NF_VALIDATE,       "VALIDATE\nTINYSA NF",menu_measure_acb},
 { MT_ADV_CALLBACK,            M_NF_AMPLIFIER,      "MEASURE\nAMP NF",menu_measure_acb},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};
#endif

static const menuitem_t menu_measure2[] = {
  { MT_ADV_CALLBACK,            M_AM,           "AM",           menu_measure_acb},
  { MT_ADV_CALLBACK,            M_FM,           "FM",           menu_measure_acb},
  { MT_ADV_CALLBACK,            M_THD,          "THD",           menu_measure_acb},
#ifdef __CHANNEL_POWER__
  { MT_ADV_CALLBACK,            M_CP,           "CHANNEL\nPOWER",menu_measure_acb},
#endif
#ifdef __LINEARITY__
{ MT_ADV_CALLBACK | MT_LOW,   M_LINEARITY,  "LINEAR",         menu_measure_acb},
#endif
#ifdef __NOISE_FIGURE__
{ MT_SUBMENU | MT_LOW,          0,            "NOISE\nFIGURE",    menu_measure_noise_figure},
#endif
#ifdef __FFT_DECONV__
  { MT_ADV_CALLBACK,            M_DECONV,  "DECONV",         menu_measure_acb},
#endif
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_measure[] = {
  { MT_ADV_CALLBACK,            M_OFF,        "OFF",            menu_measure_acb},
  { MT_ADV_CALLBACK,            M_IMD,        "HARMONIC",       menu_measure_acb},
  { MT_ADV_CALLBACK,            M_OIP3,       "OIP3",           menu_measure_acb},
  { MT_ADV_CALLBACK,            M_PHASE_NOISE,"PHASE\nNOISE",   menu_measure_acb},
  { MT_ADV_CALLBACK,            M_SNR,        "SNR",            menu_measure_acb},
  { MT_ADV_CALLBACK,            M_PASS_BAND,  "-3dB\nWIDTH",     menu_measure_acb},
  { MT_SUBMENU,  0,             S_RARROW" MORE",                menu_measure2},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

#ifdef __CALIBRATE__
#ifdef TINYSA4
static const menuitem_t menu_calibrate_harmonic[] =
{
  { MT_FORM | MT_TITLE,      0, "Connect 5.34GHz at -50 to -10dBm",  NULL},
#ifdef TINYSA4
  { MT_FORM | MT_CALLBACK,   3, "CALIBRATE",        menu_calibrate_cb},
#endif
  { MT_FORM | MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_calibrate_normal[] =
{
  { MT_FORM | MT_TITLE,      0, "Connect CAL and RF",  NULL},
#ifdef TINYSA4
  { MT_FORM | MT_CALLBACK,   1, "CALIBRATE",        menu_calibrate_cb},
#endif
  { MT_FORM | MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_calibrate[] =
{
  { MT_FORM | MT_SUBMENU,   1, "CALIBRATE 100kHz to 5.34GHz",   menu_calibrate_normal},
  { MT_FORM | MT_SUBMENU,   1, "CALIBRATE above 5.34GHz",       menu_calibrate_harmonic},
  { MT_FORM | MT_CALLBACK,   2, "RESET CALIBRATION",            menu_calibrate_cb},
  { MT_FORM | MT_NONE,     0, NULL, menu_back} // next-> menu_back
};

#else
static const menuitem_t menu_calibrate[] =
{
  { MT_FORM | MT_TITLE,      0, "Connect HIGH and LOW",  NULL},
  { MT_FORM | MT_CALLBACK,   1, "CALIBRATE",                 menu_calibrate_cb},
  { MT_FORM | MT_CALLBACK,   2, "RESET CALIBRATION",         menu_calibrate_cb},
  { MT_FORM | MT_NONE,     0, NULL, menu_back} // next-> menu_back
};
#endif

#endif

#ifdef __USE_SERIAL_CONSOLE__
const menuitem_t menu_serial_speed[] = {
  { MT_ADV_CALLBACK|MT_REPEATS, DATA_STARTS_REPEATS(0,10), "%u", menu_serial_speed_acb },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_connection[] = {
  { MT_ADV_CALLBACK, _MODE_USB,    "USB",    menu_connection_acb },
  { MT_ADV_CALLBACK, _MODE_SERIAL, "SERIAL", menu_connection_acb },
  { MT_SUBMENU,  0, "SERIAL\nSPEED", menu_serial_speed },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};
#endif

const menuitem_t menu_touch[] = {
  { MT_CALLBACK, CONFIG_MENUITEM_TOUCH_CAL,  "TOUCH CAL",  menu_config_cb},
  { MT_CALLBACK, CONFIG_MENUITEM_TOUCH_TEST, "TOUCH TEST", menu_config_cb},
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};

#ifdef __USE_RTC__
const menuitem_t menu_date_time[] = {
  { MT_KEYPAD, KM_RTC_TIME,  "SET TIME\n\b%s",                          0 },        // MUST BE BEFORE DATE
  { MT_KEYPAD, KM_RTC_DATE,  "SET DATE\n\b%s",                          0 },
  { MT_NONE, 0, NULL, menu_back} // next-> menu_back
};
#endif

static const menuitem_t menu_config2[] =
{
 { MT_ADV_CALLBACK,     0,     "PULSE\nHIGH",            menu_settings_pulse_acb},
 { MT_ADV_CALLBACK | MT_LOW, 0,"LO OUTPUT", menu_lo_output_acb},
#ifdef __ULTRA__
 { MT_ADV_CALLBACK,     0,     "ENABLE\nULTRA",    menu_ultra_acb},
#endif
 { MT_KEYPAD,   KM_GRIDLINES,  "MINIMUM\nGRIDLINES", "Enter minimum horizontal grid divisions"},
 { MT_KEYPAD,  KM_VAR,         "JOG STEP\n\b%s","0 = AUTO"},
 { MT_CALLBACK,        0 ,     "CLEAR\nCONFIG",    menu_clearconfig_cb},
#ifdef __USE_SERIAL_CONSOLE__
 { MT_SUBMENU,          0, "CONNECTION", menu_connection},
#endif
 { MT_SUBMENU,     0,              "LEVEL\nCORRECTION",    menu_actual_power},
#ifdef TINYSA4
 { MT_SUBMENU,          0, "EXPERT\nCONFIG", menu_settings},
#else
 { MT_SUBMENU,          0, "EXPERT\nCONFIG", menu_settings2},
#endif
 { MT_NONE,             0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_config[] = {
  { MT_SUBMENU,  0,                        "TOUCH",     menu_touch},
  { MT_CALLBACK, CONFIG_MENUITEM_SELFTEST, "SELF TEST", menu_config_cb},
#ifdef __CALIBRATE__
  { MT_SUBMENU,  0,                        "LEVEL CAL", menu_calibrate},
#endif
  { MT_CALLBACK, CONFIG_MENUITEM_VERSION,  "VERSION",   menu_config_cb},
#ifdef __SPUR__
  { MT_ADV_CALLBACK,0,          "%s",          menu_spur_acb},
#endif
  { MT_KEYPAD, KM_REPEAT,       "SAMPLE REP\n\b%s",    "1..100"},
#ifdef __LCD_BRIGHTNESS__
  { MT_CALLBACK, 0, "BRIGHTNESS", menu_brightness_cb},
#endif
#ifdef __USE_RTC__
  { MT_SUBMENU,  0, "DATE\nTIME", menu_date_time},
#endif
#ifndef TINYSA4
  { MT_SUBMENU,  0, S_RARROW" DFU",  menu_dfu},
#endif
  { MT_SUBMENU,  0, S_RARROW"MORE", menu_config2},
  { MT_NONE,     0, NULL, menu_back} // next-> menu_back
};
#if 0
static const menuitem_t menu_storage[] =
{
 { MT_ADV_CALLBACK,0,          "TRACE %d",        menu_storage_acb},
 { MT_ADV_CALLBACK,1,          "%s",              menu_storage_acb},
 { MT_ADV_CALLBACK,1,          "DISPLAY",         menu_storage_acb},
 { MT_ADV_CALLBACK,2,          "COPY\nFROM",      menu_storage_acb},
 { MT_ADV_CALLBACK,3,          "SUBTRACT",        menu_storage_acb},
 { MT_ADV_CALLBACK,4,          "NORMALIZE",       menu_storage_acb},
 { MT_ADV_CALLBACK,5,          "WRITE\n"S_RARROW"SD",menu_storage_acb},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};
#endif
static const menuitem_t menu_trace[] =
{
 { MT_ADV_CALLBACK|MT_REPEATS,DATA_STARTS_REPEATS(0,TRACES_MAX),          "TRACE %d",        menu_trace_acb},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_marker_trace[] =
{
 { MT_ADV_CALLBACK|MT_REPEATS,DATA_STARTS_REPEATS(0,TRACES_MAX),          "TRACE %d",        menu_marker_trace_acb},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_store_trace[] =
{
 { MT_ADV_CALLBACK|MT_REPEATS,DATA_STARTS_REPEATS(0,TRACES_MAX),          MT_CUSTOM_LABEL,        menu_store_trace_acb},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_subtract_trace[] =
{
 { MT_ADV_CALLBACK|MT_REPEATS,DATA_STARTS_REPEATS(0,TRACES_MAX+1),          MT_CUSTOM_LABEL,        menu_subtract_trace_acb},
 { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_traces[] =
{
 { MT_ADV_CALLBACK,0,          "TRACE %d",                  menu_traces_acb},
 { MT_ADV_CALLBACK,1,          "ENABLE",                    menu_traces_acb},
 { MT_ADV_CALLBACK,2,          "FREEZE",                    menu_traces_acb},
 { MT_ADV_CALLBACK,3,          MT_CUSTOM_LABEL,             menu_traces_acb},       // Calc
 { MT_ADV_CALLBACK,4,          "NORMALIZE",                 menu_traces_acb},
 { MT_ADV_CALLBACK,5,          MT_CUSTOM_LABEL,             menu_traces_acb},       // Trace Math
 { MT_SUBMENU,     0,          "COPY\n"S_RARROW"TRACE",     menu_store_trace},
#ifdef TINYSA4
 { MT_ADV_CALLBACK,6,          "WRITE\n"S_RARROW"SD",       menu_traces_acb},
#endif
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_display[] = {
  { MT_ADV_CALLBACK,0,             "PAUSE\nSWEEP",    menu_pause_acb},
  { MT_ADV_CALLBACK,1,             MT_CUSTOM_LABEL,   menu_waterfall_acb},
#ifdef __LEVEL_METER__
  { MT_ADV_CALLBACK,1,             "BIG\nNUMBER",     menu_level_meter_acb},
#endif
#ifdef __DRAW_LINE__
  { MT_ADV_CALLBACK,1,             "DRAW\nLINE",      menu_settings_draw_line_acb},
#endif
  { MT_KEYPAD,      KM_SWEEP_TIME, "SWEEP\nTIME",     "0..600s, 0=disable"},       // This must be item 3 to match highlighting
  { MT_SUBMENU,     0,             "SWEEP\nPOINTS",   menu_sweep_points},
  { MT_SUBMENU,     0,             "SWEEP\nACCURACY",  menu_sweep_speed},
  { MT_ADV_CALLBACK,0,             "ROTATE\nDISPLAY",    menu_flip_acb},
//#ifdef __REMOTE_DESKTOP__
//  { MT_ADV_CALLBACK,0,          "SEND\nDISPLAY",    menu_send_display_acb},
//#endif
//  { MT_KEYPAD,  KM_SWEEP_TIME,  "SWEEP\nTIME",    NULL},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_unit[] =
{
  { MT_ADV_CALLBACK,U_DBM,   "dBm",             menu_unit_acb},
  { MT_ADV_CALLBACK,U_DBMV,  "dBmV",            menu_unit_acb},
  { MT_ADV_CALLBACK,U_DBUV,  "dB"S_MICRO"V",    menu_unit_acb},
  { MT_ADV_CALLBACK,U_VOLT,  "Volt",            menu_unit_acb},
//{ MT_ADV_CALLBACK,U_UVOLT, S_MICRO"Volt",     menu_unit_acb},
  { MT_ADV_CALLBACK,U_WATT,  "Watt",            menu_unit_acb},
//{ MT_ADV_CALLBACK,U_UWATT, S_MICRO"Watt",     menu_unit_acb},
#ifdef TINYSA4
  { MT_ADV_CALLBACK,U_RAW,   "RAW",             menu_unit_acb},
#endif
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_trigger[] = {
  { MT_ADV_CALLBACK, T_AUTO,     "AUTO",           menu_trigger_acb},
  { MT_ADV_CALLBACK, T_NORMAL,   "NORMAL",         menu_trigger_acb},
  { MT_ADV_CALLBACK, T_SINGLE,   "SINGLE",         menu_trigger_acb},
//  { MT_ADV_CALLBACK, T_DONE,     "READY",          menu_trigger_acb},
  { MT_KEYPAD,       KM_TRIGGER, "TRIGGER LEV\n\b%s", NULL},
  { MT_ADV_CALLBACK, T_UP,       "UP\nEDGE",       menu_trigger_acb},
  { MT_ADV_CALLBACK, T_DOWN,     "DOWN\nEDGE",     menu_trigger_acb},
  { MT_ADV_CALLBACK, T_MODE,     "%s\nTRIGGER",     menu_trigger_acb},
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_level[] = {
  { MT_SUBMENU, 0,              "REF LEVEL",    menu_reflevel},
//{ MT_SUBMENU, 0,              "SCALE/DIV",    menu_scale_per},
  { MT_CALLBACK,0,              "SCALE/DIV",    menu_scale_cb},
  { MT_SUBMENU, 0,              "ATTENUATE",    menu_atten},
//  { MT_SUBMENU,0,             "CALC",         menu_average},
  { MT_SUBMENU, 0,              "UNIT",         menu_unit},
  { MT_KEYPAD,  KM_EXT_GAIN,      "EXT GAIN\n\b%s",NULL},
#ifdef TINYSA4
  { MT_ADV_CALLBACK | MT_LOW ,0,"LNA",          menu_extra_lna_acb},
 #endif
  { MT_SUBMENU,  0,             "TRIGGER",      menu_trigger},
#ifdef __LISTEN__
  { MT_ADV_CALLBACK, 0,             "LISTEN",       menu_listen_acb},
#endif
  { MT_NONE,   0, NULL, menu_back} // next-> menu_back
};

static const menuitem_t menu_stimulus[] = {
  { MT_KEYPAD,  KM_START,       "START\n\b%s",       NULL},
  { MT_KEYPAD,  KM_STOP,        "STOP\n\b%s",        NULL},
  { MT_KEYPAD,  KM_CENTER,      "CENTER\n\b%s",      NULL},
  { MT_KEYPAD,  KM_SPAN,        "SPAN\n\b%s",        NULL},
  { MT_KEYPAD,  KM_CW,          "ZERO SPAN",   NULL},
  { MT_SUBMENU,0,               "RBW",         menu_rbw},
#ifdef __VBW__
  { MT_SUBMENU,     0,             "VBW",             menu_vbw},
#endif
  { MT_ADV_CALLBACK,0,          "SHIFT\nFREQ", menu_shift_acb},
  { MT_NONE,    0, NULL, menu_back} // next-> menu_back
};

#ifdef TINYSA4
const menuitem_t menu_mode[] = {
//  { MT_FORM | MT_TITLE,                 0,                      "tinySA MODE",           NULL},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_LOW_INPUT+I_SA,       "Spectrum Analyzer",      menu_mode_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_LOW_OUTPUT+I_SINUS,   "Signal Generator",     menu_mode_acb},
//  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_HIGH_OUTPUT+I_GEN,    "%s to HIGH out",    menu_mode_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_CONNECT+I_GEN,        "Calibration Output: %s",   menu_sreffer_acb},
  { MT_FORM | MT_NONE,     0, NULL, NULL } // sentinel
};
#else
const menuitem_t menu_mode[] = {
//  { MT_FORM | MT_TITLE,                 0,                      "tinySA MODE",           NULL},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_LOW_INPUT+I_SA,       "%s to LOW in",      menu_mode_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_HIGH_INPUT+I_SA,      "%s to HIGH in",     menu_mode_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_LOW_OUTPUT+I_SINUS,   "%s to LOW out",     menu_mode_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_HIGH_OUTPUT+I_GEN,    "%s to HIGH out",    menu_mode_acb},
  { MT_FORM | MT_ADV_CALLBACK | MT_ICON,    I_CONNECT+I_GEN,        "Cal. output: %s",   menu_sreffer_acb},
//  { MT_SUBMENU,  0, "EXPERT\nCONFIG", menu_settings3},
//  { MT_FORM | MT_CANCEL,   0, S_RARROW" BACK", NULL },
  { MT_FORM | MT_NONE,     0, NULL, NULL } // sentinel
};
#endif

static const menuitem_t menu_top[] = {
  { MT_SUBMENU,  0, "PRESET",       menu_load_preset},
  { MT_SUBMENU,  0, "FREQUENCY",    menu_stimulus},
  { MT_SUBMENU,  0, "LEVEL",        menu_level},
  { MT_SUBMENU,  0, "TRACE",        menu_traces},
  { MT_SUBMENU,  0, "DISPLAY",      menu_display},
  { MT_SUBMENU,  0, "MARKER",       menu_marker},
  { MT_SUBMENU,  0, "MEASURE",      menu_measure},
  { MT_SUBMENU,  0, "CONFIG",       menu_config},
  { MT_SUBMENU,  0, "MODE",         menu_mode},
  { MT_NONE,     0, NULL, NULL } // sentinel,
 // MENUITEM_CLOSE,
};

// ===[MENU DEFINITION END]======================================================

#define ACTIVE_COLOR RGBHEX(0x007FFF)

static void menu_item_modify_attribute(                     // To modify menu buttons with keypad modes
    const menuitem_t *menu, int item, ui_button_t *button)
{
  if (menu == menu_display) {
    if (item == 5)
      button->icon = setting.sweep_time_us != 0 ? BUTTON_ICON_CHECK_MANUAL : BUTTON_ICON_CHECK_AUTO;
//  } else if (menu == menu_sweep_speed) {
//    if (item == 3)
//    button->icon = setting.fast_speedup != 0 ? BUTTON_ICON_CHECK_MANUAL : BUTTON_ICON_CHECK_AUTO;
  } else if (menu == menu_reflevel) {
    if (item == 1)
      button->icon = setting.auto_reflevel ? BUTTON_ICON_GROUP: BUTTON_ICON_GROUP_CHECKED;
  } else if (menu == menu_atten) {
    if (item == 1)
      button->icon = setting.auto_attenuation ? BUTTON_ICON_GROUP: BUTTON_ICON_GROUP_CHECKED;
  }
}

static void fetch_numeric_target(uint8_t mode)
{
  switch (mode) {
  case KM_START:
    uistat.freq_value = get_sweep_frequency(ST_START) + (setting.frequency_offset - FREQUENCY_SHIFT);
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_STOP:
    uistat.freq_value = get_sweep_frequency(ST_STOP) + (setting.frequency_offset - FREQUENCY_SHIFT);
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_CENTER:
    uistat.freq_value = get_sweep_frequency(ST_CENTER) + (setting.frequency_offset - FREQUENCY_SHIFT);
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_SPAN:
    uistat.freq_value = get_sweep_frequency(ST_SPAN);
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_CW:
    uistat.freq_value = get_sweep_frequency(ST_CW) + (setting.frequency_offset - FREQUENCY_SHIFT);
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_SCALE:
  case KM_LINEAR_SCALE:
    uistat.value = setting.scale;
    plot_printf(uistat.text, sizeof uistat.text, "%f/", uistat.value);
    break;
  case KM_REFLEVEL:
    uistat.value = setting.reflevel;
    plot_printf(uistat.text, sizeof uistat.text, "%.3F", uistat.value);
    break;
  case KM_ATTENUATION:
    uistat.value = get_attenuation();
    plot_printf(uistat.text, sizeof uistat.text, "%ddB", ((int32_t)uistat.value));
     break;
  case KM_ACTUALPOWER:
    uistat.value = get_level_offset();
    plot_printf(uistat.text, sizeof uistat.text, "%ddB", ((int32_t)uistat.value));
    break;
  case KM_IF:
    uistat.freq_value = setting.frequency_IF;
    if (!setting.auto_IF)
      plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    else
      plot_printf(uistat.text, sizeof uistat.text, "AUTO");
    break;
#ifdef TINYSA4
  case KM_IF2:
    uistat.freq_value = config.frequency_IF2;
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_R:
    uistat.value = setting.R;
    if (setting.R)
      plot_printf(uistat.text, sizeof uistat.text, "%d", setting.R);
    else
      plot_printf(uistat.text, sizeof uistat.text, "AUTO");
    break;
  case KM_MOD:
    if (local_modulo)
      plot_printf(uistat.text, sizeof uistat.text, "%d", local_modulo);
    else
      plot_printf(uistat.text, sizeof uistat.text, "AUTO");
   break;
  case KM_CP:
    uistat.value = ADF4350_modulo;
    plot_printf(uistat.text, sizeof uistat.text, "%d", ADF4350_modulo);
    break;
#endif
  case KM_SAMPLETIME:
    uistat.value = setting.step_delay;
    if (uistat.value)
      plot_printf(uistat.text, sizeof uistat.text, "%dus", ((int32_t)uistat.value));
    else
      plot_printf(uistat.text, sizeof uistat.text, "AUTO");
    break;
  case KM_FAST_SPEEDUP:
    uistat.value = setting.fast_speedup;
    plot_printf(uistat.text, sizeof uistat.text, "%d", ((int32_t)uistat.value));
    break;
  case KM_REPEAT:
    uistat.value = setting.repeat;
    plot_printf(uistat.text, sizeof uistat.text, "%d", ((int32_t)uistat.value));
    break;
  case KM_LOWOUTLEVEL:
    uistat.value = get_level();           // compensation for dB offset during low output mode
    float end_level =  ((int32_t)uistat.value)+setting.level_sweep;
    if (end_level < level_min())
      end_level = level_min();
    if (end_level > level_max())
      end_level = level_max();
    uistat.value += setting.external_gain;
    end_level += setting.external_gain;
    if (setting.level_sweep != 0)
      plot_printf(uistat.text, sizeof uistat.text, "%.1f to %.1fdBm", uistat.value, end_level);
    else
#ifdef TINYSA4
      plot_printf(uistat.text, sizeof uistat.text, "%+.1fdBm %s", uistat.value, (setting.disable_correction?"Uncorrected":""));
#else
      plot_printf(uistat.text, sizeof uistat.text, "%+.1fdBm", uistat.value);
#endif
    break;
  case KM_HIGHOUTLEVEL:
    uistat.value = get_level();           // compensation for dB offset during low output mode
    uistat.value += setting.external_gain;
    plot_printf(uistat.text, sizeof uistat.text, "%+.1fdBm", uistat.value);
    break;
  case KM_DECAY:
    uistat.value = setting.decay;
    plot_printf(uistat.text, sizeof uistat.text, "%d", ((int32_t)uistat.value));
    break;
#ifdef __QUASI_PEAK__
  case KM_ATTACK:
    uistat.value = setting.attack;
    plot_printf(uistat.text, sizeof uistat.text, "%d", ((int32_t)uistat.value));
    break;
#endif
#ifdef __ULTRA__
  case KM_ULTRA_START:
    uistat.freq_value = config.ultra_start;
    if (config.ultra_start == ULTRA_AUTO)
      plot_printf(uistat.text, sizeof uistat.text, "AUTO");
    else
      plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value );
    break;
  case KM_DIRECT_START:
    uistat.freq_value = config.direct_start;
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_DIRECT_STOP:
    uistat.freq_value = config.direct_stop;
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
#endif
#ifdef __LIMITS__
  case KM_LIMIT_FREQ:
    uistat.freq_value = setting.limits[current_trace][active_limit].frequency;
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_LIMIT_LEVEL:
    uistat.value = value(setting.limits[current_trace][active_limit].level);
    plot_printf(uistat.text, sizeof uistat.text, "%.1f", uistat.value);
    break;
#endif
  case KM_NOISE:
    uistat.value = setting.noise;
    plot_printf(uistat.text, sizeof uistat.text, "%d", ((int32_t)uistat.value));
    break;
#ifdef TINYSA4
  case KM_FREQ_CORR:
    if (config.setting_frequency_30mhz >= 3000000000ULL)
      uistat.value = (config.setting_frequency_30mhz - 3000000000ULL)/3;
    else
      uistat.value = - ((int)(3000000000ULL - config.setting_frequency_30mhz))/3;

    plot_printf(uistat.text, sizeof uistat.text, "%d", (int32_t)uistat.value);
    break;
#else
  case KM_10MHZ:
    uistat.freq_value = config.setting_frequency_10mhz;
    plot_printf(uistat.text, sizeof uistat.text, "%3.6fMHz", uistat.freq_value);
    break;
#endif
  case KM_EXT_GAIN:
    uistat.value = setting.external_gain;
    plot_printf(uistat.text, sizeof uistat.text, "%.1fdB", uistat.value);
    break;
  case KM_LEVELSWEEP:
    uistat.value = setting.level_sweep;
    plot_printf(uistat.text, sizeof uistat.text, "%.1fdB", uistat.value);
    break;
  case KM_SWEEP_TIME:
//    if (setting.sweep_time_us < calc_min_sweep_time_us())
//      uistat.value = calc_min_sweep_time_us();
//    else
      uistat.value = setting.sweep_time_us;
    uistat.value /= (float)ONE_SECOND_TIME;
    plot_printf(uistat.text, sizeof uistat.text, "%.3Fs", uistat.value);
    break;
  case KM_TRIGGER:
    uistat.value = value(setting.trigger_level);
    char *format;
    if (UNIT_IS_LINEAR(setting.unit))
      format = "%.3F%s"; // 5 characters incl u, m, etc...
    else
      format = "%.1f%s";
    plot_printf(uistat.text, sizeof uistat.text, format, uistat.value,unit_string[setting.unit]);
    break;
  case KM_MARKER:
    if (active_marker >=0) {
      uistat.freq_value = markers[active_marker].frequency;
      plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    }
    break;
  case KM_MODULATION:
    if (active_marker >=0) {
      uistat.value = setting.modulation_frequency;
      plot_printf(uistat.text, sizeof uistat.text, "%7.0fHz", uistat.value);
    }
    break;
#ifdef TINYSA4
  case KM_DEVIATION:
    uistat.freq_value = setting.modulation_deviation_div100 * 100;
    plot_printf(uistat.text, sizeof uistat.text, "%.3QHz", uistat.freq_value);
    break;
  case KM_DEPTH:
    uistat.value = setting.modulation_depth_x100;
    plot_printf(uistat.text, sizeof uistat.text, "%3d", (int)uistat.value);
    break;
#endif
  case KM_VAR:
    uistat.freq_value = setting.frequency_var;
    if ( setting.frequency_var)
      plot_printf(uistat.text, sizeof uistat.text, "%.4QHz", setting.frequency_var);
    else
      plot_printf(uistat.text, sizeof uistat.text, "AUTO");

    break;
#ifdef __NOISE_FIGURE__
  case KM_NF:
    uistat.value = config.noise_figure;
    plot_printf(uistat.text, sizeof uistat.text, "%.1fdB", uistat.value);
    break;
#endif
#ifdef __USE_RTC__
  case KM_RTC_TIME:
  {
    uint32_t tr = rtc_get_tr_bin(); // TR read first
    plot_printf(uistat.text, sizeof uistat.text, "%02d:%02d:%02d",
      RTC_TR_HOUR(dr),
      RTC_TR_MIN(dr),
      RTC_TR_SEC(dr));
  }
  break;
  case KM_RTC_DATE:
  {
    uint32_t dr = rtc_get_dr_bin(); // DR read second
    plot_printf(uistat.text, sizeof uistat.text, "20%02d/%02d/%02d",
      RTC_DR_YEAR(dr),
      RTC_DR_MONTH(dr),
      RTC_DR_DAY(dr));
  }
    break;
#endif
  }
}

static void
set_numeric_value(void)
{
  switch (keypad_mode) {
  case KM_START:
    set_sweep_frequency(ST_START, uistat.freq_value - (setting.frequency_offset - FREQUENCY_SHIFT));
    break;
  case KM_STOP:
    set_sweep_frequency(ST_STOP, (freq_t)(uistat.freq_value - (setting.frequency_offset - FREQUENCY_SHIFT)));
    break;
  case KM_CENTER:
    set_sweep_frequency(ST_CENTER, uistat.freq_value - (setting.frequency_offset - FREQUENCY_SHIFT));
    break;
  case KM_SPAN:
    setting.modulation = MO_NONE;
    set_sweep_frequency(ST_SPAN, uistat.freq_value);
    break;
  case KM_CW:
    set_sweep_frequency(ST_CW, uistat.freq_value - (setting.frequency_offset - FREQUENCY_SHIFT));
    break;
  case KM_LINEAR_SCALE:
  case KM_SCALE:
    user_set_scale(uistat.value);
    break;
  case KM_REFLEVEL:
    user_set_reflevel(uistat.value);
    break;
  case KM_ATTENUATION:
    setting.auto_attenuation = false;
    set_attenuation(uistat.value);
    break;
  case KM_ACTUALPOWER:
    set_actual_power(uistat.value);
    config_save();
    break;
  case KM_IF:
    set_IF(uistat.freq_value);
//    config_save();
    break;
#ifdef TINYSA4
  case KM_IF2:
    set_IF2(uistat.freq_value);
//    config_save();
    break;
  case KM_R:
    set_R(uistat.value);
//    config_save();
    break;
  case KM_MOD:
    set_modulo(uistat.value);
    break;
  case KM_CP:
    ADF4351_CP((int)uistat.value);
//    config_save();
    break;
#endif
  case KM_SAMPLETIME:
    set_step_delay(uistat.value);
    break;
  case KM_OFFSET_DELAY:
    set_offset_delay(uistat.value);
    break;
  case KM_FAST_SPEEDUP:
    set_fast_speedup(uistat.value);
    break;
  case KM_REPEAT:
    set_repeat(uistat.value);
    break;
  case KM_LOWOUTLEVEL:
    set_level(uistat.value - setting.external_gain);
    break;
  case KM_HIGHOUTLEVEL:
    set_level(uistat.value - setting.external_gain);
    break;
  case KM_DECAY:
    set_decay(uistat.value);
    break;
#ifdef __QUASI_PEAK__
  case KM_ATTACK:
    set_attack(uistat.value);
    break;
#endif
#ifdef __ULTRA__
  case KM_ULTRA_START:
    config.ultra_start = uistat.freq_value;
    reset_settings(setting.mode);
//    config_save(); // TODO not now
    //ultra_start = config.ultra_start;
    break;
  case KM_DIRECT_START:
    config.direct_start = uistat.freq_value;
    config_save();
    break;
  case KM_DIRECT_STOP:
    config.direct_stop = uistat.freq_value;
    config_save();
    break;
#endif
#ifdef TINYSA4
  case KM_EXP_AVER:
    setting.exp_aver = uistat.value;
    dirty = true;
#endif
  case KM_LEVEL:
    break;
#ifdef __LIMITS__
  case KM_LIMIT_FREQ:
    setting.limits[current_trace][active_limit].frequency = uistat.freq_value - (setting.frequency_offset - FREQUENCY_SHIFT);
    dirty = true;
    limits_update();
    break;
  case KM_LIMIT_LEVEL:
    setting.limits[current_trace][active_limit].level = to_dBm(uistat.value);
    dirty = true;
    limits_update();
    break;
#endif
  case KM_NOISE:
    set_noise(uistat.value);
    break;
#ifdef TINYSA4
  case KM_FREQ_CORR:
    set_freq_corr(uistat.value);
    break;
#else
  case KM_10MHZ:
    set_10mhz(uistat.freq_value);
    break;
#endif
  case KM_EXT_GAIN:
    set_external_gain(uistat.value);
    break;
  case KM_LEVELSWEEP:
    setting.modulation = MO_NONE;
    set_level_sweep(uistat.value);
    break;
  case KM_SWEEP_TIME:
    set_sweep_time_us(uistat.value*ONE_SECOND_TIME);
    update_grid();
    break;
  case KM_TRIGGER:
    if (setting.trigger == T_AUTO )
      set_trigger(T_NORMAL);
    set_trigger_level(to_dBm(uistat.value));
    completed = true;
    break;
  case KM_GRIDLINES:
    set_gridlines(uistat.value);
    break;
  case KM_MARKER:
    set_marker_frequency(active_marker, uistat.freq_value - (setting.frequency_offset - FREQUENCY_SHIFT));
    break;
  case KM_MARKER_TIME:
    set_marker_time(active_marker, uistat.value);
    break;
  case KM_MODULATION:
    set_modulation_frequency((int)uistat.value);
    break;
#ifdef TINYSA4
  case KM_COR_AM:
    config.cor_am = -(int)uistat.value;
    config_save();
    dirty = true;
    break;
  case KM_COR_WFM:
    config.cor_wfm = -(int)uistat.value;
    config_save();
    dirty = true;
    break;
  case KM_COR_NFM:
    config.cor_nfm = -(int)uistat.value;
    config_save();
    dirty = true;
    break;
  case KM_DEVIATION:
    set_deviation((int)uistat.freq_value);
    break;
  case KM_DEPTH:
    set_depth((int)uistat.value);
    break;
#endif
  case KM_VAR:
    setting.frequency_var = uistat.freq_value;
    break;
#ifdef __NOISE_FIGURE__
  case KM_NF:
    config.noise_figure = uistat.value;
    config_save();
    dirty = true;
    break;
#endif
#ifdef __USE_RTC__
  case KM_RTC_DATE:
  case KM_RTC_TIME:
    {
      int i = 0;
      uint32_t  dt_buf[2];
      dt_buf[0] = rtc_get_tr_bcd(); // TR should be read first for sync
      dt_buf[1] = rtc_get_dr_bcd(); // DR should be read second
      //            0    1   2       4      5     6
      // time[] ={sec, min, hr, 0, day, month, year, 0}
      uint8_t   *time = (uint8_t*)dt_buf;
      for (; i < 6 && kp_buf[i]!=0; i++) kp_buf[i]-= '0';
      for (; i < 6                ; i++) kp_buf[i] =   0;
      for (i = 0; i < 3; i++) kp_buf[i] = (kp_buf[2*i]<<4) | kp_buf[2*i+1]; // BCD format
      if (keypad_mode == KM_RTC_DATE) {
        // Month limit 1 - 12 (in BCD)
             if (kp_buf[1] <    1) kp_buf[1] =    1;
        else if (kp_buf[1] > 0x12) kp_buf[1] = 0x12;
        // Day limit (depend from month):
        uint8_t day_max = 28 + ((0b11101100000000000010111110111011001100>>(kp_buf[1]<<1))&3);
        day_max = ((day_max/10)<<4)|(day_max%10); // to BCD
             if (kp_buf[2] <  1)      kp_buf[2] = 1;
        else if (kp_buf[2] > day_max) kp_buf[2] = day_max;
        time[6] = kp_buf[0]; // year
        time[5] = kp_buf[1]; // month
        time[4] = kp_buf[2]; // day
      }
      else {
        // Hour limit 0 - 23, min limit 0 - 59, sec limit 0 - 59 (in BCD)
        if (kp_buf[0] > 0x23) kp_buf[0] = 0x23;
        if (kp_buf[1] > 0x59) kp_buf[1] = 0x59;
        if (kp_buf[2] > 0x59) kp_buf[2] = 0x59;
        time[2] = kp_buf[0]; // hour
        time[1] = kp_buf[1]; // min
        time[0] = kp_buf[2]; // sec
      }
      rtc_set_time(dt_buf[1], dt_buf[0]);
    }
    break;
#endif

  }
}

void
menu_move_top(void)
{
  while (menu_current_level > 0)
    menu_move_back(false);
}


// -------------------------- CAL STATUS ---------------------------------------------
const char * const dBText[] = { "1dB/", "2dB/", "5dB/", "10dB/", "20dB/"};
const int refMHz[] = { 30, 15, 10, 4, 3, 2, 1 };

float my_round(float v)
{
  float m = 1;
  int sign = 1;
  if (v < 0) {
    sign = -1;
    v = -v;
  }
  while (v < 100) {
    v = v * 10;
    m = m / 10;
  }
  while (v > 1000) {
    v = v / 10;
    m = m * 10;
  }
  v = (int)(v+0.5);
  v = v * m;
  if (sign == -1) {
    v = -v;
  }
  return v;
}
const char * const unit_string[MAX_UNIT_TYPE*2] = { "dBm", "dBmV", "dB"S_MICRO"V", "RAW", "V", "W", "dB", "dB", "dB", "RAW", "V", "W" }; // unit + 6 is delta unit

//static const float scale_value[]={50000, 20000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20,10,5,2,1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002, 0.001,0.0005,0.0002, 0.0001};
//static const char * const scale_vtext[]= {"50000", "20000", "10000", "5000", "2000", "1000", "500", "200", "100", "50", "20","10","5","2","1","0.5","0.2","0.1","0.05","0.02","0.01", "0.005","0.002","0.001", "0.0005","0.0002","0.0001"};

// Quick menu
#define MAX_QUICK_MENU  20
#define MAX_ITEM_SPACE   2
static uint16_t    quick_menu_y[MAX_QUICK_MENU];
static menuitem_t  *quick_menu[MAX_QUICK_MENU];
static uint8_t max_quick_menu = 0;
static uint8_t item_space = 0; //

int invoke_quick_menu(int y)
{
  int i;
  for (i = 0; i < max_quick_menu;i++) {
    if (y < quick_menu_y[i]) {
      if ((uint32_t)quick_menu[i] < KM_NONE) {
        ui_mode_keypad((int)quick_menu[i]);
      } else {
        selection = -1;
        menu_current_level = 0;
        menu_push_submenu(quick_menu[i]);
      }
      return TRUE;
    }
  }
  return FALSE;
}
#define YSTEP   8

int add_quick_menu(int y, menuitem_t *menu)
{
  y += YSTEP*item_space/2 + YSTEP;
  if (max_quick_menu<MAX_QUICK_MENU-1) {
    quick_menu_y[max_quick_menu] = y;
    quick_menu[max_quick_menu++] = menu;
  }
  return y;
}

const char *month[] = { "Jan", "Feb", "Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec" };

void draw_cal_status(void)
{
#define BLEN    7
  char buf[BLEN+1];
  buf[6]=0;
  int x = 0;
  int y = OFFSETY;
  unsigned int color;
  const bool rounding = !UNIT_IS_LINEAR(setting.unit);
  const char * const unit = unit_string[setting.unit];
redraw_cal_status:
  buf[6]=0;
  x = 0;
  y = OFFSETY;
  ili9341_set_background(LCD_BG_COLOR);
  ili9341_fill(0, 0, OFFSETX, LCD_HEIGHT);
  max_quick_menu = 0;
  if (MODE_OUTPUT(setting.mode)) {     // No cal status during output
#ifdef TINYSA4
    if (level_error) {
      ili9341_set_foreground(LCD_BRIGHT_COLOR_RED);
      ili9341_drawstring("LEVEL\nERROR", 0 , 80);
    }
    if (depth_error) {
      ili9341_set_foreground(LCD_BRIGHT_COLOR_RED);
      ili9341_drawstring("DEPTH\nERROR", 0 , 140);
    }
#endif
    return;
  }

    //  if (current_menu_is_form() && !in_selftest)
//    return;

//ili9341_set_background(LCD_BG_COLOR);

  float yMax = setting.reflevel;

  if (level_is_calibrated())
    color = setting.auto_reflevel ? LCD_FG_COLOR : LCD_BRIGHT_COLOR_GREEN;
  else
    color = LCD_BRIGHT_COLOR_RED;
  ili9341_set_foreground(color);
  // Ref level
  if (rounding)
    lcd_printf(x, y, "%+4d", (int)yMax);
  else
    lcd_printf(x, y, "%+4.3F", (yMax/setting.unit_scale));
  y = add_quick_menu(y, (menuitem_t *)menu_reflevel);

  // Unit
#if 0
  color = LCD_FG_COLOR;
  ili9341_set_foreground(color);
  if (setting.auto_reflevel){
    y += YSTEP + YSTEP/2 ;
    ili9341_drawstring("AUTO", x, y);
  }
#endif
  lcd_printf(x, y, "%c%s", unit_scale_text[setting.unit_scale_index], unit);
  y = add_quick_menu(y, (menuitem_t *)menu_unit);

  // Scale
  ili9341_set_foreground(LCD_FG_COLOR);
#if 0
  unsigned int i = 0;
  while (i < ARRAY_COUNT(scale_value)) {
    float t = (setting.scale/setting.unit_scale) / scale_value[i];
    if (t > 0.9 && t < 1.1){
      lcd_printf(x, y, "%s%c/",scale_vtext[i],unit_scale_text[setting.unit_scale_index]);
      break;
    }
    i++;
  }
#else
  lcd_printf(x, y, "%.2F/",setting.scale);
#endif
  y = add_quick_menu(y, (menuitem_t *)KM_SCALE);

  // Trigger status
  if (is_paused()) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    ili9341_drawstring("PAUSED", x, y);
    y += YSTEP + YSTEP/2 ;
  }
  if (setting.trigger == T_SINGLE || setting.trigger == T_NORMAL ) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    ili9341_drawstring("ARMED", x, y);
    y += YSTEP + YSTEP/2 ;
  }
  // AM warning
  if (signal_is_AM) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_RED);
    ili9341_drawstring("AM", x, y);
    y += YSTEP + YSTEP/2 ;
  }
  quick_menu_y[max_quick_menu] = y;
  quick_menu[max_quick_menu++] = (menuitem_t *)NULL;

//  if (setting.mode == M_LOW) {
    // Attenuation
    ili9341_set_foreground(setting.auto_attenuation ? LCD_FG_COLOR : LCD_BRIGHT_COLOR_GREEN);
    lcd_printf(x, y, "Atten:\n%4.2FdB", get_attenuation());
    y = add_quick_menu(y+= YSTEP, (menuitem_t *)menu_atten);
//  }

  // Calc
  if (setting.average[0]>0) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    lcd_printf(x, y, "Calc:\n%s", averageText[setting.average[0]]);
    y = add_quick_menu(y+= YSTEP, (menuitem_t *)menu_average);
  }
#ifdef __SPUR__  // Spur
#ifdef TINYSA3
  if (setting.spur_removal != S_OFF) {
#endif
    ili9341_set_foreground(setting.spur_removal == S_ON ? LCD_BRIGHT_COLOR_GREEN : LCD_FG_COLOR);
    lcd_printf(x, y, "Spur:\n%s", S_IS_AUTO(setting.spur_removal) ? "AUTO" : (setting.spur_removal == S_OFF ?"OFF" : "ON"));
    y = add_quick_menu(y += YSTEP, (menuitem_t *)menu_config);
#ifdef TINYSA3
  }
#endif
  if (setting.mirror_masking) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    ili9341_drawstring("Mask:\nON", x, y);
    y = add_quick_menu(y+=YSTEP, (menuitem_t *)menu_stimulus);
  }
#endif

  if (setting.subtract[0]) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    ili9341_drawstring("Norm.", x, y);
    y = add_quick_menu(y, (menuitem_t *)menu_display);
  }

  // RBW
  ili9341_set_foreground(setting.rbw_x10 ? LCD_BRIGHT_COLOR_GREEN : LCD_FG_COLOR);
  if (dirty) update_rbw();
  lcd_printf(x, y, "RBW:\n%.1FHz", actual_rbw_x10*100.0);
  y = add_quick_menu(y+=YSTEP, (menuitem_t *)menu_rbw);

#ifdef __VBW__
  // VBW
  if (setting.frequency_step > 0) {
    int vbw = setting.vbw_x100;
    if (vbw != 0)
      color = LCD_BRIGHT_COLOR_GREEN;
    else {
      color = LCD_FG_COLOR;
      vbw = 1;
    }
    ili9341_set_foreground(color);
    lcd_printf(x, y, "VBW:\n%.1FHz", actual_rbw_x10*100.0 / vbw);
    y = add_quick_menu(y+=YSTEP, (menuitem_t *)menu_vbw);
  }
#endif
  // Sweep time: SD_NORMAL, SD_PRECISE, SD_FAST, SD_MANUAL
  static const char fscan[]={0, 'P', 'F', 'N', 'M'};
  if (dirty) {
    calculate_step_delay();
    setting.actual_sweep_time_us = calc_min_sweep_time_us();
  }
#if 0                   // Activate for sweep time debugging
  lcd_printf(x, y, "%cScan:\n%5.3Fs", fscan[setting.step_delay_mode&7], (float)setting.sweep_time_us/ONE_SECOND_TIME);
#endif
  ili9341_set_foreground((setting.step_delay_mode&7) != 0 ? LCD_BRIGHT_COLOR_GREEN : LCD_FG_COLOR);
  lcd_printf(x, y, "%cScan:", fscan[setting.step_delay_mode&7]);
  ili9341_set_foreground((setting.step_delay || setting.sweep_time_us ) ? LCD_BRIGHT_COLOR_GREEN : LCD_FG_COLOR);
  lcd_printf(x, y+YSTEP, "%5.3Fs",(float)setting.actual_sweep_time_us/ONE_SECOND_TIME);
  y = add_quick_menu(y+=YSTEP, (menuitem_t *)menu_sweep_speed);


  #if 0                   // Activate for sweep time debugging
  y += YSTEP;
  update_rbw();             // To ensure the calc_min_sweep time shown takes the latest delay into account
  calculate_step_delay();
  uint32_t t = calc_min_sweep_time_us();
  lcd_printf(x, y, "%5.3Fs", (float)t/ONE_SECOND_TIME);
  y += YSTEP;
  lcd_printf(x, y, "%5.3Fs", (float)setting.additional_step_delay_us/ONE_SECOND_TIME);
  y += YSTEP + YSTEP/2 ;
#endif
#ifdef TINYSA4
  if (setting.extra_lna){
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    lcd_printf(x, y, "LNA:ON");
    y = add_quick_menu(y, (menuitem_t *)menu_level);
    y += YSTEP;
  }

  if (config.ultra_start != ULTRA_AUTO){
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    lcd_printf(x, y, "Ultra:\n%3QHz", ultra_start);
    y = add_quick_menu(y += YSTEP, (menuitem_t *)menu_config);
  }

  if (setting.disable_correction){
    ili9341_set_foreground(LCD_BRIGHT_COLOR_RED);
    lcd_printf(x, y, "Corr:\nOFF");
    y += 2*YSTEP + YSTEP/2;
  }

  if (force_signal_path){
    ili9341_set_foreground(LCD_BRIGHT_COLOR_RED);
    lcd_printf(x, y, "Path:\n%s", path_text[signal_path]);
    y += 2*YSTEP + YSTEP/2;
  }

#endif
  // Cal output
  if (setting.refer >= 0) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    lcd_printf(x, y, "Ref:\n%dMHz",reffer_freq[setting.refer]/1000000);
    y = add_quick_menu(y+=YSTEP, (menuitem_t *)menu_reffer);
  }

  // Offset
  if (setting.external_gain != 0.0) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    lcd_printf(x, y, "Gain:\n%4.1fdB",setting.external_gain);
    y = add_quick_menu(y+=YSTEP, (menuitem_t *)KM_EXT_GAIN);
  }

  // Repeat
  if (setting.repeat != 1) {
    ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    lcd_printf(x, y, "Repeat\n x%d", setting.repeat);
    y = add_quick_menu(y+=YSTEP,( menuitem_t *)KM_REPEAT);
  }

  // Trigger
  if (setting.trigger != T_AUTO) {
    if (is_paused() || setting.trigger == T_NORMAL) {
      ili9341_set_foreground(LCD_BRIGHT_COLOR_GREEN);
    } else {
      ili9341_set_foreground(LCD_BRIGHT_COLOR_RED);
    }
    ili9341_drawstring("TRIG:", x, y);

    y += YSTEP;
    if (rounding)
      lcd_printf(x, y, "%6.3f", value(setting.trigger_level));
    else
      lcd_printf(x, y, "%6.4F", value(setting.trigger_level));
//    lcd_printf(x, y, "%4f", value(setting.trigger_level)/setting.unit_scale);
    y = add_quick_menu(y,(menuitem_t *)menu_trigger);
  }
#ifndef TINYSA4
  // Mode
  ili9341_set_foreground(level_is_calibrated() ? LCD_BRIGHT_COLOR_GREEN : LCD_BRIGHT_COLOR_RED);
  ili9341_drawstring_7x13(MODE_LOW(setting.mode) ? "LOW" : "HIGH", x, y);
  y += YSTEP + YSTEP/2 ;
#endif
  // Compact status string
//  ili9341_set_background(LCD_FG_COLOR);
  ili9341_set_foreground(LCD_FG_COLOR);
  strncpy(buf,"     ",BLEN-1);
  if (setting.auto_IF)
    buf[0] = 'f';
  else
    buf[0] = 'F';
  if (S_IS_AUTO(setting.agc))
    buf[1] = 'g';
  else if (S_STATE(setting.agc))
    buf[1] = 'G';
  if (S_IS_AUTO(setting.lna))
    buf[2] = 'n';
  else if (S_STATE(setting.lna))
    buf[2] = 'N';
  if (S_IS_AUTO(setting.below_IF))
    buf[3] = 'b';
  else if (S_STATE(setting.below_IF))
    buf[3] = 'B';
#ifdef TINYSA4
  if (S_IS_AUTO(setting.spur_removal))
    buf[4] = 's';
  else if (S_STATE(setting.spur_removal))
    buf[4] = 'S';
#endif
  ili9341_drawstring(buf, x, y);

  // Version
  y += YSTEP + YSTEP/2 ;
#ifdef TINYSA4 // 'tinySA4_v1.2-[0-9]*-gxxxxxxx'
  strncpy(buf,&TINYSA_VERSION[9], BLEN+1); // '1.2-...'
#else // 'tinySA_v1.2-[0-9]*-gxxxxxxx'
  strncpy(buf,&TINYSA_VERSION[8], BLEN+1); // '1.2-...'
#endif
  if (buf[5]=='-' ) { // '1.2-n-g...'
    if (buf[4]=='0')  // '1.2-0-g...'
      buf[3] = 0;  // -> '1.2'
    else {
      buf[5] = buf[4]; // -> '1.200n'
      buf[4] = '0';
      buf[3] = '0';
    }
  } else if (buf[6]=='-' ) { // 1.2-nn-g...
    buf[3] = '0'; // -> '1.20nn'
  } else { // 1.2-345-g... (or 1.2-3456...)
    buf[3] = buf[4]; // -> '1.2345'
    buf[4] = buf[5];
    buf[5] = buf[6];
  }
  buf[6] = 0;
  ili9341_drawstring(buf, x, y);

#ifdef __USE_RTC__
  y += YSTEP + YSTEP/2 ;
  uint32_t dr = rtc_get_dr_bin(); // DR read second
  lcd_printf(x, y,  "20%02d/\n%s/%02d", RTC_DR_YEAR(dr), month[RTC_DR_MONTH(dr)-1], RTC_DR_DAY(dr));
  y += YSTEP*2;
  uint32_t tr = rtc_get_tr_bin(); // TR read first
  lcd_printf(x, y,  "%02d:%02d", RTC_TR_HOUR(dr), RTC_TR_MIN(dr));
#endif


  if (y >= BATTERY_START && item_space > 0) {
    item_space--;                       // Reduce item spacing
    goto redraw_cal_status;
  }
  if ((y + (max_quick_menu+1) * YSTEP/2) < BATTERY_START && item_space < MAX_ITEM_SPACE) {
    item_space++;                       // Increase item spacing
    goto redraw_cal_status;
  }

//  ili9341_set_background(LCD_BG_COLOR);
  if (!setting.waterfall) {               // Do not draw bottom level if in waterfall mode
    // Bottom level
    y = area_height + OFFSETY;
    if (level_is_calibrated())
      if (setting.auto_reflevel)
        color = LCD_FG_COLOR;
      else
        color = LCD_BRIGHT_COLOR_GREEN;
    else
      color = LCD_BRIGHT_COLOR_RED;
    ili9341_set_foreground(color);
    if (rounding)
      lcd_printf(x, y, "%4d", (int)(yMax - setting.scale * NGRIDY));
    else
      lcd_printf(x, y, "%+4.3F", ((yMax - setting.scale * NGRIDY)/setting.unit_scale));
    y = add_quick_menu(y,(menuitem_t *)menu_average);
  }
}

//----------------------------

#define MENU_STACK_DEPTH_MAX 7
const menuitem_t *menu_stack[MENU_STACK_DEPTH_MAX] = {
  menu_top, NULL, NULL, NULL, NULL, NULL, NULL
};

int current_menu_is_form(void)
{
  return menu_stack[menu_current_level]->type & MT_FORM;
}

static bool menuDisabled(uint8_t type){
  if ((type & MT_LOW) && !MODE_LOW(setting.mode))
    return true;
  if ((type & MT_HIGH) && !MODE_HIGH(setting.mode))
    return true;
//  if (type == MT_BLANK)
//    return true;
  return false;
}

static const menuitem_t *menu_next_item(const menuitem_t *m, int *sub_item){
  do{
    if (m->type & MT_REPEATS) {
      (*sub_item)++;
      if (*sub_item < ((m->data>>4) & 0x0f))
        return m;
      *sub_item = 0;
    }
    m++;
    m = MT_MASK(m->type) == MT_NONE ? (menuitem_t *)m->reference : m;
  } while(m!=NULL && menuDisabled(m->type));
  return m;
}

static const menuitem_t *current_menu_item(int i, int *sub_item){
  *sub_item = 0;
  const menuitem_t * m = menu_stack[menu_current_level];
  while (i--) m = menu_next_item(m,sub_item);
  return m;
}

static int current_menu_get_count(void){
  int i = 0,sub_item = 0;
  const menuitem_t *m = menu_stack[menu_current_level];
  while (m){m = menu_next_item(m, &sub_item); i++;}
  return i;
}

static void
ensure_selection(void)
{
  const menuitem_t *menu = menu_stack[menu_current_level];
  int i = current_menu_get_count();
  if (selection <  0) selection =  -1;
  if (selection >= i) selection = i-1;
  if (MT_MASK(menu[0].type) == MT_TITLE && selection == 0) selection = 1;
  if (i <  MENU_BUTTON_MIN) i = MENU_BUTTON_MIN;
  if (i >= MENU_BUTTON_MAX) i = MENU_BUTTON_MAX;
#ifdef MENU_USE_AUTOHEIGHT
  menu_button_height = MENU_BUTTON_HEIGHT_N(i);
#endif
}

static void
menu_move_back(bool leave_ui)
{
  if (menu_current_level == 0)
    return;
  erase_menu_buttons();
  bool form = current_menu_is_form();
  menu_current_level--;

  // redraw all if switch from form to normal menu mode or back
  if (form != current_menu_is_form())
    redraw_request|=REDRAW_AREA|REDRAW_BATTERY|REDRAW_FREQUENCY|REDRAW_CAL_STATUS;

  selection = -1;
  if (leave_ui)
    ui_mode_normal();
  else
    ui_mode_menu();
}

void
menu_push_submenu(const menuitem_t *submenu)
{
  erase_menu_buttons();
  if (menu_current_level < MENU_STACK_DEPTH_MAX-1)
    menu_current_level++;
  menu_stack[menu_current_level] = submenu;
  ui_mode_menu();
}

void
menu_push_lowoutput(void)
{
  menu_push_submenu(menu_lowoutputmode);
}

void
menu_push_highoutput(void)
{
  menu_push_submenu(menu_highoutputmode);
}

/*
static void
menu_move_top(void)
{
  if (menu_current_level == 0)
    return;
  menu_current_level = 0;
  ensure_selection();
  erase_menu_buttons();
  draw_menu();
}
*/

static void
menu_invoke(int item)
{
  int sub_item;
  const menuitem_t *menu = current_menu_item(item, &sub_item);
  if (menu == NULL) return;
  switch (MT_MASK(menu->type)) {
//  case MT_NONE:
//  case MT_BLANK:
//    ui_mode_normal();
//    break;

  case MT_CANCEL:
    menu_move_back(false);
    break;

  case MT_CALLBACK: {
    uistat.auto_center_marker = false;
    menuaction_cb_t cb = (menuaction_cb_t)menu->reference;
    if (cb) (*cb)(item, (menu->type & MT_REPEATS) ? (menu->data & 0x0f)+sub_item : menu->data);
//    if (!(menu->type & MT_FORM))
    redraw_request |= REDRAW_CAL_STATUS;
    break;
  }
  case MT_ADV_CALLBACK: {
    uistat.auto_center_marker = false;
    menuaction_acb_t cb = (menuaction_acb_t)menu->reference;
    if (cb) (*cb)(item, (menu->type & MT_REPEATS) ? (menu->data & 0x0f)+sub_item : menu->data, NULL);
//    if (!(menu->type & MT_FORM))
    redraw_request |= REDRAW_CAL_STATUS | REDRAW_BATTERY;
    break;
  }
  case MT_SUBMENU:
    menu_push_submenu((const menuitem_t*)menu->reference);
    break;

  case MT_KEYPAD:
    uistat.auto_center_marker = false;
    if (current_menu_is_form()) {
      redraw_frame();         // Remove form numbers
    }
    kp_help_text = (char *)menu->reference;
    if (menu->data <= KM_CW) {      // One of the frequency input keypads
      if (MODE_LOW(setting.mode))
        kp_help_text = VARIANT("0..350MHz",range_text);
      else if (menu->data == KM_SPAN)
        kp_help_text = VARIANT("0..720Mhz",range_text);
      else
        kp_help_text = VARIANT("240..960Mhz",range_text);
    }
    ui_mode_keypad(menu->data);
    redraw_request |= REDRAW_CAL_STATUS;
    break;
  }
  // Redraw menu after if UI in menu mode
  if (ui_mode == UI_MENU)
    draw_menu();
}

static const char * const keypad_scale_text[] = {"0", "1", "2", "5", "10", "20" , "50", "100", "200", "500"};
//static const int  keypad_scale_value[] = { 1, 2, 5, 10, 20 , 50, 100, 200, 500};

static void
draw_button(uint16_t x, uint16_t y, uint16_t w, uint16_t h, ui_button_t *b)
{
  uint16_t bw = b->border&BUTTON_BORDER_WIDTH_MASK;
  ili9341_set_foreground(b->fg);
  ili9341_set_background(b->bg);ili9341_fill(x + bw, y + bw, w - (bw * 2), h - (bw * 2));
  if (bw==0) return;
  uint16_t br = LCD_RISE_EDGE_COLOR;
  uint16_t bd = LCD_FALLEN_EDGE_COLOR;
  uint16_t type = b->border;
  ili9341_set_background(type&BUTTON_BORDER_TOP    ? br : bd);ili9341_fill(x,          y,           w, bw); // top
  ili9341_set_background(type&BUTTON_BORDER_RIGHT  ? br : bd);ili9341_fill(x + w - bw, y,          bw,  h); // right
  ili9341_set_background(type&BUTTON_BORDER_LEFT   ? br : bd);ili9341_fill(x,          y,          bw,  h); // left
  ili9341_set_background(type&BUTTON_BORDER_BOTTOM ? br : bd);ili9341_fill(x,          y + h - bw,  w, bw); // bottom
  // Set colors for button text after
  ili9341_set_background(b->bg);
}

void drawMessageBox(char *header, char *text, uint32_t delay){
  ui_button_t b;
  b.bg = LCD_MENU_COLOR;
  b.fg = LCD_MENU_TEXT_COLOR;
  b.border = BUTTON_BORDER_FLAT|1;
  // Draw header
  draw_button((LCD_WIDTH-MESSAGE_BOX_WIDTH)/2, LCD_HEIGHT/2-40, MESSAGE_BOX_WIDTH, 60, &b);
  ili9341_drawstring_7x13(header, (LCD_WIDTH-MESSAGE_BOX_WIDTH)/2 + 10, LCD_HEIGHT/2-40 + 5);
  // Draw window
  ili9341_set_background(LCD_FG_COLOR);
  ili9341_fill((LCD_WIDTH-MESSAGE_BOX_WIDTH)/2+3, LCD_HEIGHT/2-40+bFONT_STR_HEIGHT+8, MESSAGE_BOX_WIDTH-6, 60-bFONT_STR_HEIGHT-8-3);
  ili9341_drawstring_7x13(text, (LCD_WIDTH-MESSAGE_BOX_WIDTH)/2 + 20, LCD_HEIGHT/2-40 + bFONT_STR_HEIGHT + 8 + 14);
  chThdSleepMilliseconds(delay);
}

static void
draw_keypad(uint32_t mask)
{
  int i;
  ui_button_t button;
  button.fg = LCD_MENU_TEXT_COLOR;
  for(i = 0; keypads[i].c >= 0; i++) {
    if ((mask&(1<<i)) == 0) continue;
    button.bg = LCD_MENU_COLOR;
    if (i == selection){
      button.bg = LCD_MENU_ACTIVE_COLOR;
      button.border = KEYBOARD_BUTTON_BORDER|BUTTON_BORDER_FALLING;
    }
    else
      button.border = KEYBOARD_BUTTON_BORDER|BUTTON_BORDER_RISE;
    int x = KP_GET_X(keypads[i].x);
    int y = KP_GET_Y(keypads[i].y);
    draw_button(x, y, KP_WIDTH, KP_HEIGHT, &button);
    if (keypads[i].c < KP_0) { // KP_0
      ili9341_drawfont(keypads[i].c,
                     x + (KP_WIDTH - NUM_FONT_GET_WIDTH) / 2,
                     y + (KP_HEIGHT - NUM_FONT_GET_HEIGHT) / 2);
    } else {
      const char *t = keypad_scale_text[keypads[i].c - KP_0];
      ili9341_drawstring_10x14(t,
                     x + (KP_WIDTH  - wFONT_MAX_WIDTH*strlen(t)) / 2,
                     y + (KP_HEIGHT - wFONT_GET_HEIGHT) / 2);
    }
  }
}

static int
menu_is_multiline(const char *label)
{
  int n = 1;
  if (label)
    while (*label)
      if (*label++ == '\n')
        n++;
  return n;
}

static int period_pos(void) {int j; for (j = 0; j < kp_index && kp_buf[j] != '.'; j++); return j;}

static void
draw_numeric_input(const char *buf)
{
  uint16_t i;
  uint16_t x = 10 + 10 * FONT_WIDTH + 4;
  uint16_t xsim;
#ifdef __USE_RTC__
  if (keypad_mode == KM_RTC_DATE || keypad_mode == KM_RTC_TIME)
    xsim = 0b01010100;
  else
#endif
  xsim = (0b00100100100100100 >>(2-(period_pos()%3)))&(~1);
  ili9341_set_foreground(LCD_INPUT_TEXT_COLOR);
  ili9341_set_background(LCD_INPUT_BG_COLOR);
  for (i = 0; buf[i]; i++) {
    int c = buf[i];
         if (c == '.'){c = KP_PERIOD;xsim<<=4;}
    else if (c == '-'){c = KP_MINUS; xsim&=~3;}
    else// if (c >= '0' && c <= '9')
      c = c - '0';
    if (c < 0) c = 0;
    // Add space before char
    int16_t space = xsim&1 ? 2 + 10 : 2;
    xsim>>=1;
    ili9341_fill(x, LCD_HEIGHT-NUM_INPUT_HEIGHT+4, space, NUM_INPUT_HEIGHT);
    x+=space;
    if (c >= 0) // c is number
      ili9341_drawfont(c, x, LCD_HEIGHT-NUM_INPUT_HEIGHT+4);
    else
      break;
    x+=NUM_FONT_GET_WIDTH;
  }

  ili9341_fill(x, LCD_HEIGHT-NUM_INPUT_HEIGHT+4, LCD_WIDTH - 1 - x, NUM_INPUT_HEIGHT);
  if (buf[0] == 0 && kp_help_text != NULL) {
    int  lines = menu_is_multiline(kp_help_text);
    ili9341_set_foreground(LCD_INPUT_TEXT_COLOR);
    ili9341_drawstring_7x13(kp_help_text, 64+NUM_FONT_GET_WIDTH+2, LCD_HEIGHT-(lines*bFONT_GET_HEIGHT+NUM_INPUT_HEIGHT)/2);
  }
}

static void
draw_numeric_area_frame(void)
{
  ili9341_set_foreground(LCD_INPUT_TEXT_COLOR);
  ili9341_set_background(LCD_INPUT_BG_COLOR);

  ili9341_fill(0, LCD_HEIGHT-NUM_INPUT_HEIGHT, LCD_WIDTH, NUM_INPUT_HEIGHT);
  char *name = keypads_mode_tbl[keypad_mode].name;
  int lines = menu_is_multiline(name);
  ili9341_drawstring_7x13(name, 10, LCD_HEIGHT-NUM_INPUT_HEIGHT + (NUM_INPUT_HEIGHT-lines*bFONT_STR_HEIGHT)/2);
  //ili9341_drawfont(KP_KEYPAD, 300, 216);
  draw_numeric_input("");
}

#define ICON_WIDTH        16
#define ICON_HEIGHT       11

static const uint8_t check_box[] = {
  _BMP16(0b0011111111110000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0011111111110000),

  _BMP16(0b0011111111110000),
  _BMP16(0b0010000000001000),
  _BMP16(0b0010000000011000),
  _BMP16(0b0010000000110000),
  _BMP16(0b0010000001100000),
  _BMP16(0b0010100011010000),
  _BMP16(0b0010110110010000),
  _BMP16(0b0010011100010000),
  _BMP16(0b0010001000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0011111111110000),

  _BMP16(0b0011111111111000),
  _BMP16(0b0010000000001000),
  _BMP16(0b0010001111101000),
  _BMP16(0b0010011001101000),
  _BMP16(0b0010110001101000),
  _BMP16(0b0010110001101000),
  _BMP16(0b0010111111101000),
  _BMP16(0b0010110001101000),
  _BMP16(0b0010110001101000),
  _BMP16(0b0010000000001000),
  _BMP16(0b0011111111111000),

  _BMP16(0b0011111111111000),
  _BMP16(0b0010000000001000),
  _BMP16(0b0010110001101000),
  _BMP16(0b0010110001101000),
  _BMP16(0b0010111011101000),
  _BMP16(0b0010111111101000),
  _BMP16(0b0010110101101000),
  _BMP16(0b0010110101101000),
  _BMP16(0b0010110001101000),
  _BMP16(0b0010000000001000),
  _BMP16(0b0011111111111000),

  _BMP16(0b0000000000000000),
  _BMP16(0b0000011110000000),
  _BMP16(0b0000100001000000),
  _BMP16(0b0001000000100000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0010000000010000),
  _BMP16(0b0001000000100000),
  _BMP16(0b0000100001000000),
  _BMP16(0b0000011110000000),

  _BMP16(0b0000000000000000),
  _BMP16(0b0000011110000000),
  _BMP16(0b0000100001000000),
  _BMP16(0b0001001100100000),
  _BMP16(0b0010011110010000),
  _BMP16(0b0010111111010000),
  _BMP16(0b0010111111010000),
  _BMP16(0b0010011110010000),
  _BMP16(0b0001001100100000),
  _BMP16(0b0000100001000000),
  _BMP16(0b0000011110000000),
};

#ifndef MENU_USE_AUTOHEIGHT
#ifdef TINYSA4
#define menu_button_height  ((menu[i].type & MT_FORM) || menu_is_multiline(menu[i].label) == 2 ? LCD_HEIGHT/10 : LCD_HEIGHT/12 )
#else
#define menu_button_height  ((menu[i].type & MT_FORM) || menu_is_multiline(menu[i].label) == 2 ? LCD_HEIGHT/8 : LCD_HEIGHT/10 )
#endif
#endif

static void
draw_menu_buttons(const menuitem_t *menu, uint32_t mask)
{
  int i, y;
  ui_button_t button;
  const menuitem_t *m = menu;
  int sub_item = 0;
//  while (menuDisabled(m->type))
//    m = menu_next_item(m, &sub_item);       // Just in case the first item is disabled
  for (i = 0, y = 0; m; m = menu_next_item(m, &sub_item), i++, y += menu_button_height) {
    if ((mask&(1<<i)) == 0)
      continue;
    button.icon = BUTTON_ICON_NONE;
    // Border width
    button.border = MENU_BUTTON_BORDER;

    if (MT_MASK(m->type) == MT_TITLE) {
      button.fg = LCD_FG_COLOR;
      button.bg = LCD_BG_COLOR;
      button.border = 0; // no border for title
    } else {
      button.bg = LCD_MENU_COLOR;
      button.fg = LCD_MENU_TEXT_COLOR;
    }

    if (i == selection){
      button.bg = LCD_MENU_ACTIVE_COLOR;
      button.border|= BUTTON_BORDER_FALLING;
    }
    else
      button.border|= BUTTON_BORDER_RISE;

    // Need replace this obsolete bad function on new MT_ADV_CALLBACK variant
    menu_item_modify_attribute(menu, i, &button);      // before plot_printf to create status text
    char *text;
    // MT_ADV_CALLBACK - allow change button data in callback, more easy and correct
    if (MT_MASK(m->type) == MT_ADV_CALLBACK){
      menuaction_acb_t cb = (menuaction_acb_t)m->reference;
      if (cb) (*cb)(i, (m->type & MT_REPEATS) ? (m->data & 0x0f)+sub_item : m->data, &button);
      // Apply custom text, from button label and
      if (m->label != MT_CUSTOM_LABEL)
        plot_printf(button.text, sizeof(button.text), m->label, button.param_1.u);
      text = button.text;
    }
    else
      text = (char *)m->label;
    // Only keypad retrieves value
    if (MT_MASK(m->type) == MT_KEYPAD) {
      fetch_numeric_target(m->data);
      plot_printf(button.text, sizeof button.text, m->label, uistat.text);
      text = button.text;
    }

    int button_height = menu_button_height;
    if (current_menu_is_form()) {
      int button_width = MENU_FORM_WIDTH;
      int button_start = (LCD_WIDTH - MENU_FORM_WIDTH)/2; // At center of screen
      draw_button(button_start, y, button_width, button_height, &button);
      uint16_t text_offs = button_start + 6;
      if (button.icon >=0){
        ili9341_blitBitmap(button_start+3, y+(button_height-ICON_HEIGHT)/2, ICON_WIDTH, ICON_HEIGHT, &check_box[button.icon*2*ICON_HEIGHT]);
        text_offs = button_start+6+ICON_WIDTH+1;
      }
#ifdef __ICONS__
      if (m->type & MT_ICON) {
        ili9341_blitBitmap(button_start+MENU_FORM_WIDTH-2*FORM_ICON_WIDTH-8,y+(button_height-FORM_ICON_HEIGHT)/2,FORM_ICON_WIDTH,FORM_ICON_HEIGHT,& left_icons[((menu[i].data >>4)&0xf)*2*FORM_ICON_HEIGHT]);
        ili9341_blitBitmap(button_start+MENU_FORM_WIDTH-  FORM_ICON_WIDTH-8,y+(button_height-FORM_ICON_HEIGHT)/2,FORM_ICON_WIDTH,FORM_ICON_HEIGHT,&right_icons[((menu[i].data >>0)&0xf)*2*FORM_ICON_HEIGHT]);
      }
#endif
      int local_text_shift = 0;
      if (MT_MASK(m->type) == MT_KEYPAD) {
        int local_slider_positions = 0;
        if (m->data == KM_CENTER) {
          local_slider_positions =  LCD_WIDTH/2+setting.slider_position;
          lcd_printf(button_start+12 + 0 * MENU_FORM_WIDTH/5, y+button_height-9, "%+3.0FHz", -(float)setting.slider_span);
          lcd_printf(button_start+12 + 1 * MENU_FORM_WIDTH/5, y+button_height-9, "%+3.0FHz", -(float)setting.slider_span/10);
          lcd_printf(button_start+12 + 2 * MENU_FORM_WIDTH/5, y+button_height-9, "Set");
          lcd_printf(button_start+12 + 3 * MENU_FORM_WIDTH/5, y+button_height-9, "%+3.0FHz",  (float)setting.slider_span/10);
          lcd_printf(button_start+12 + 4 * MENU_FORM_WIDTH/5, y+button_height-9, "%+3.0FHz",  (float)setting.slider_span);
        } else if (m->data == KM_LOWOUTLEVEL) {
          local_slider_positions = ((get_level() - level_min()) * (MENU_FORM_WIDTH-8)) / level_range() + OFFSETX+4;
          lcd_printf(button_start+12 + 0 * MENU_FORM_WIDTH/5, y+button_height-9, "%+ddB", -10);
          lcd_printf(button_start+12 + 1 * MENU_FORM_WIDTH/5, y+button_height-9, "%+ddB",  -1);
          lcd_printf(button_start+12 + 2 * MENU_FORM_WIDTH/5, y+button_height-9, "Set");
          lcd_printf(button_start+12 + 3 * MENU_FORM_WIDTH/5, y+button_height-9, "%+ddB",   1);
          lcd_printf(button_start+12 + 4 * MENU_FORM_WIDTH/5, y+button_height-9, "%+ddB",  10);
        } else if (m->data == KM_HIGHOUTLEVEL) {
          local_slider_positions = ((get_level() - level_min() ) * (MENU_FORM_WIDTH-8)) / level_range() + OFFSETX+4;
        }
        if (local_slider_positions){
          if (local_slider_positions < button_start)
            local_slider_positions = button_start;
          // Shift down text if slider present
          local_text_shift = 2;
          ili9341_blitBitmap(local_slider_positions - 4, y, 7, 5, slider_bitmap);
          // Draw divider for sliders
          for (int i = 1; i <= 4; i++)
            ili9341_line(button_start + i * MENU_FORM_WIDTH/5, y+button_height-9, button_start + i * MENU_FORM_WIDTH/5, y+button_height);
        }
      }
//      ili9341_drawstring_size(text, text_offs, y+(button_height-2*FONT_GET_HEIGHT)/2-local_text_shift, 2);
      ili9341_drawstring_10x14(text, text_offs, y+(button_height-wFONT_GET_HEIGHT)/2-local_text_shift);
    } else {
      int button_width = MENU_BUTTON_WIDTH;
      int button_start = LCD_WIDTH - MENU_BUTTON_WIDTH;
      draw_button(button_start, y, button_width, button_height, &button);
      uint16_t text_offs = button_start + 7;
      if (button.icon >=0){
        ili9341_blitBitmap(button_start+2, y+(button_height-ICON_HEIGHT)/2, ICON_WIDTH, ICON_HEIGHT, &check_box[button.icon*2*ICON_HEIGHT]);
        text_offs = button_start+2+ICON_WIDTH;
      }
      int lines = menu_is_multiline(text);
#define BIG_BUTTON_FONT 1
#ifdef BIG_BUTTON_FONT
      ili9341_drawstring_7x13(text, text_offs, y+(button_height-lines*bFONT_GET_HEIGHT)/2);
#else
      ili9341_drawstring(text, text_offs, y+(button_height-linesFONT_GET_HEIGHT)/2);
#endif
    }
  }
  // Cleanup other buttons (less flicker)
  // Erase empty buttons
  if (NO_WATERFALL - y > 0){
    ili9341_set_background(LCD_BG_COLOR);
    ili9341_fill(LCD_WIDTH-MENU_BUTTON_WIDTH, y, MENU_BUTTON_WIDTH, NO_WATERFALL - y);
  }
//  if (current_menu_is_form())
//    draw_battery_status();
}

enum { SL_UNKNOWN, SL_SPAN, SL_MOVE};

void set_keypad_value(int v)
{
  keypad_mode = v;
  set_numeric_value();
}

void check_frequency_slider(freq_t slider_freq)
{

  if ( (maxFreq - minFreq) < (freq_t)setting.slider_span) {
    setting.slider_span = maxFreq - minFreq;                         // absolute mode with max step size
  }
  freq_t half_span = setting.slider_span >> 1;
  int temp = (setting.slider_span / (MENU_FORM_WIDTH-8));
  if (minFreq + half_span > slider_freq) {
    setting.slider_position -= (minFreq + half_span - slider_freq) / temp;            // reposition if needed
  }
  if (maxFreq < slider_freq + half_span) {
    setting.slider_position += (slider_freq + half_span - maxFreq) / temp;            // reposition if needed
  }
}

#define TOUCH_DEAD_ZONE 40
static void
menu_select_touch(const menuitem_t * m, int i, int pos)
{
  uint32_t mask = (1<<i)|(1<<selection);
  selection = i;
  draw_menu_mask(mask);
#if 1               // drag values
  int keypad = m->data;
  int touch_x, touch_y,  prev_touch_x = 0;
  systime_t dt = 0;
  int mode = SL_UNKNOWN;

  systime_t ticks = chVTGetSystemTimeX();
  while (touch_check() != EVT_TOUCH_NONE){
    dt = chVTGetSystemTimeX() - ticks;
    if (dt > BUTTON_DOWN_LONG_TICKS) break;
  }

  if (current_menu_is_form() && MT_MASK(m->type) == MT_KEYPAD && dt >= BUTTON_DOWN_LONG_TICKS){
    // Wait release touch and process it
    while (touch_check() != EVT_TOUCH_NONE){
      touch_position(&touch_x, &touch_y);
      if (abs(touch_x -  prev_touch_x) < 2) continue;

      fetch_numeric_target(keypad);
      int new_slider = touch_x - LCD_WIDTH/2;   // Can have negative outcome
      if (new_slider < - (MENU_FORM_WIDTH-8)/2 - 1)
        new_slider = -(MENU_FORM_WIDTH-8)/2 - 1;
      if (new_slider > (MENU_FORM_WIDTH-8)/2 + 1)
        new_slider = (MENU_FORM_WIDTH-8)/2 + 1;
      if (keypad == KM_CENTER) {
        if (mode == SL_UNKNOWN ) {
          if (abs(setting.slider_position - new_slider) < TOUCH_DEAD_ZONE) // Pick up slider
            mode = SL_MOVE;
          else
            mode = SL_SPAN;
        }
        if (mode == SL_MOVE ) {
          long_t freq_delta = (setting.slider_span/(MENU_FORM_WIDTH-8))*(new_slider - setting.slider_position);
          if (freq_delta < 0 && uistat.freq_value < (freq_t)(-freq_delta))
            uistat.freq_value = 0;
          else
            uistat.freq_value+= freq_delta;
          if (uistat.freq_value < minFreq)
            uistat.freq_value = minFreq;
          if (uistat.freq_value > maxFreq)
            uistat.freq_value = maxFreq;
          setting.slider_position = new_slider;
          set_keypad_value(keypad);
          dirty = false;
          perform(false, 0, uistat.freq_value, false);
          draw_menu_mask(1<<i);
        }
        else if (mode == SL_SPAN ){
          freq_t slider_freq;
          slider_freq = uistat.freq_value;
          int pw=new_slider + LCD_WIDTH/2;
          setting.slider_position = pw - LCD_WIDTH/2;   // Show delta on slider
          setting.slider_span = 10;
          while (pw>0) {
            setting.slider_span += setting.slider_span;
            pw -= 12;
            if (pw <=0)
              break;
            setting.slider_span += setting.slider_span + (setting.slider_span >>1);
            pw -= 12;
            if (pw<=0)
              break;
            setting.slider_span *= 2;
            pw -= 12;
          }
          if ((freq_t)setting.slider_span > (maxFreq - minFreq))
            setting.slider_span = (maxFreq - minFreq);
          freq_t old_minFreq = minFreq;           // Save when in high mode
          minFreq = 0;                              // And set minFreq to 0 for span display
          uistat.freq_value = setting.slider_span;
          set_keypad_value(keypad);
          plot_printf(center_text, sizeof center_text, "RANGE: %%s");
          draw_menu_mask(1<<i);                      // Show slider span
          minFreq = old_minFreq;                     // and restore minFreq
          uistat.freq_value = slider_freq;        // and restore current slider freq
          set_keypad_value(keypad);
          plot_printf(center_text, sizeof center_text, "FREQ: %%s");
          setting.slider_position = 0;                         // reset slider after span change
          check_frequency_slider(slider_freq);
        }
        chThdSleepMilliseconds(100);
      } else if (keypad == KM_LOWOUTLEVEL) {
        uistat.value =  setting.external_gain + ((touch_x - OFFSETX+4) * level_range() ) / (MENU_FORM_WIDTH-8) + level_min() ;
        set_keypad_value(keypad);
        draw_menu_mask(1<<i);
        perform(false, 0, get_sweep_frequency(ST_CENTER), false);
      } else if (keypad == KM_HIGHOUTLEVEL) {
        set_level( (touch_x - OFFSETX+4) *(level_range()) / (MENU_FORM_WIDTH-8) + level_min() );
        draw_menu_mask(1<<i);
        perform(false, 0, get_sweep_frequency(ST_CENTER), false);
      }
      prev_touch_x = touch_x;

    }
    selection = -1;
    draw_menu_mask(1<<i);
    return;
  }
  if (current_menu_is_form() && MT_MASK(m->type) == MT_KEYPAD){
    bool do_exit = false;
    long_t step = 0;
    if (keypad == KM_LOWOUTLEVEL) {
      switch (pos) {
        case 0:step = -10;break;
        case 1:step =  -1;break;
        case 2:goto nogo;
        case 3:step =  +1;break;
        case 4:step = +10;break;
      }
      uistat.value = setting.external_gain + get_level() + step;
      do_exit = true;
    }
    else if (keypad == KM_CENTER) {
      switch (pos) {
        case 0: step = setting.slider_span;   step=-step; break;
        case 1: step = setting.slider_span/10;step=-step; break;
        case 2: goto nogo;
        case 3: step = setting.slider_span/10;break;
        case 4: step = setting.slider_span;   break;
      }
      if (step < 0 && get_sweep_frequency(ST_CENTER) < (freq_t)(-step))
        uistat.freq_value = 0;
      else
        uistat.freq_value = get_sweep_frequency(ST_CENTER) + step;
      do_exit = true;
      setting.slider_position = 0;                         // reset slider after step
      check_frequency_slider(uistat.freq_value);
    }
    if (do_exit){
      set_keypad_value(keypad);
      selection = -1;
      draw_menu_mask(1<<i);
      perform(false, 0, get_sweep_frequency(ST_CENTER), false);
      return;
    }
  }
nogo:
  setting.slider_position = 0;            // Reset slider when entering frequency
#endif
//  touch_wait_release();
  selection = -1;
  menu_invoke(i);
}


static void
menu_apply_touch(int touch_x, int touch_y)
{
  const menuitem_t *m = menu_stack[menu_current_level];
  int i;
  int y = 0;
  int active_button_start;
  if (current_menu_is_form()) {
    active_button_start = (LCD_WIDTH - MENU_FORM_WIDTH)/2;
//  active_button_stop = LCD_WIDTH - active_button_start;
  } else {
    active_button_start = LCD_WIDTH - MENU_BUTTON_WIDTH;
//  active_button_stop = LCD_WIDTH;
  }
  int sub_item = 0;
//  while (menuDisabled(m->type))
//    m = menu_next_item(m, &sub_item);       // Just in case the first item is disabled
  for (i = 0; m; m = menu_next_item(m,&sub_item), i++, y+= menu_button_height) {
    if (MT_MASK(m->type) == MT_TITLE) continue;
    if (y < touch_y && touch_y < y+menu_button_height && touch_x > active_button_start) {
      menu_select_touch(m, i, (( touch_x - active_button_start) * 5 ) / MENU_FORM_WIDTH);
      return;
    }
  }
  if (current_menu_is_form())
    return;
  touch_wait_release();
  ui_mode_normal();
}

void
draw_menu(void)
{
  draw_menu_buttons(menu_stack[menu_current_level], -1);
}

void
draw_menu_mask(uint32_t mask)
{
  draw_menu_buttons(menu_stack[menu_current_level], mask);
}

#ifdef __SWEEP_RESTART__
systime_t old_sweep_time;

void
refresh_sweep_menu(int i)
{
  current_index = i;
  systime_t new_sweep_time = chVTGetSystemTimeX();
  if (new_sweep_time - old_sweep_time > ONE_SECOND_TIME/200 && i >= 0) {
    old_sweep_time = new_sweep_time;
    if (menu_stack[menu_current_level] == menu_lowoutputmode)
      draw_menu_buttons(menu_stack[menu_current_level], 1 << 5);
    if (menu_stack[menu_current_level] == menu_highoutputmode)
      draw_menu_buttons(menu_stack[menu_current_level], 1 << 5);
  }
}
#endif

static void
erase_menu_buttons(void)
{
// Not need, screen redraw in all cases
//  ili9341_fill(area_width, 0, LCD_WIDTH - area_width, area_height, LCD_BG_COLOR);
 // if (current_menu_is_form())
 //   ili9341_fill(OFFSETX, 0,LCD_WIDTH-OFFSETX, menu_button_height*MENU_BUTTON_MAX, LCD_BG_COLOR);
 // else
 //   ili9341_fill(LCD_WIDTH-MENU_BUTTON_WIDTH, 0, MENU_BUTTON_WIDTH, menu_button_height*MENU_BUTTON_MAX, LCD_BG_COLOR);
  draw_frequencies();
}

#if 0
static void
erase_numeric_input(void)
{
  ili9341_set_background(LCD_BG_COLOR);
  ili9341_fill(0, LCD_HEIGHT-NUM_INPUT_HEIGHT, LCD_WIDTH, NUM_INPUT_HEIGHT);
}
#endif

static void
leave_ui_mode()
{
//  if (ui_mode == UI_MENU) {
//    request_to_draw_cells_behind_menu();
//    erase_menu_buttons();
//  }
  ili9341_set_background(LCD_BG_COLOR);
  // Erase bottom area (not redraw on area update)
//  if (menu_button_height*MENU_BUTTON_MAX - area_height > 0)
//    ili9341_fill(LCD_WIDTH-MENU_BUTTON_WIDTH, area_height, MENU_BUTTON_WIDTH, menu_button_height*MENU_BUTTON_MAX - area_height);
  if (setting.waterfall)
    set_waterfall();
  redraw_request|=REDRAW_AREA | REDRAW_FREQUENCY | REDRAW_CAL_STATUS | REDRAW_BATTERY;
}

void
ui_mode_menu(void)
{
//  if (ui_mode == UI_MENU)
//    return;
  ui_mode = UI_MENU;
  ensure_selection();
  if (current_menu_is_form()) {
    redraw_frame();
    area_width = 0;
    area_height = 0;
  } else {
    area_width = AREA_WIDTH_NORMAL - MENU_BUTTON_WIDTH;
    area_height = AREA_HEIGHT_NORMAL;
  }
  draw_menu();
  redraw_request|=REDRAW_BATTERY|REDRAW_CAL_STATUS;
}

static void
ui_mode_keypad(int _keypad_mode)
{
  if (ui_mode == UI_KEYPAD && keypad_mode == _keypad_mode )
    return;

  // keypads array
  keypad_mode = _keypad_mode;
  keypads = keypads_mode_tbl[_keypad_mode].keypad_type;

  ui_mode = UI_KEYPAD;
  if (!current_menu_is_form())
    draw_menu();
  draw_keypad(-1);
  draw_numeric_area_frame();
  ui_process_keypad();
}

void
ui_mode_normal(void)
{
  if (ui_mode == UI_NORMAL)
    return;
  if (current_menu_is_form())
    return;
  leave_ui_mode();
  area_width  = AREA_WIDTH_NORMAL;
  area_height = AREA_HEIGHT_NORMAL;
  ui_mode = UI_NORMAL;
}

static void
lever_move_marker(int status)
{
  uint16_t step = 1<<2;
  do {
    if (active_marker != MARKER_INVALID && markers[active_marker].enabled) {
      int idx = (int)markers[active_marker].index;
      if (status & EVT_DOWN) {
        idx -= step>>2;
        if (idx < 0) idx = 0 ;
      }
      if (status & EVT_UP) {
        idx += step>>2;
        if (idx  > sweep_points-1) idx = sweep_points-1 ;
      }
      markers[active_marker].index = idx;
      markers[active_marker].frequency = getFrequency(idx);
      redraw_marker(active_marker);
      markers[active_marker].mtype &= ~M_TRACKING;    // Disable tracking when dragging marker
      step++;
    }
    status = btn_wait_release();
  } while (status != 0);
}

static void
lever_search_marker(int status)
{
  int i = -1;
  if (active_marker != MARKER_INVALID) {
    if (status & EVT_DOWN)
      i = marker_search_left_max(markers[active_marker].index);
    else if (status & EVT_UP)
      i = marker_search_right_max(markers[active_marker].index);
    if (i != -1) {
      markers[active_marker].index = i;
      interpolate_maximum(active_marker);
//      markers[active_marker].frequency = frequencies[i];
    }
    redraw_marker(active_marker);
  }
}

// ex. 10942 -> 10000
//      6791 ->  5000
//       341 ->   200
static freq_t
step_round(freq_t v)
{
  // decade step
  freq_t x = 1;
  for (x = 1; x*10 <= v; x*= 10)
    ;

  // 1-2-5 step
  if (x * 2 > v)
    return x;
  else if (x * 5 > v)
    return x * 2;
  else
    return x * 5;
}

static void
lever_zoom_time(int status)
{
  uint32_t time = setting.sweep_time_us; // in uS
  if (time < MINIMUM_SWEEP_TIME)
    time = MINIMUM_SWEEP_TIME;
  if (status & EVT_UP) {
    time = time*10/25;
  } else if (status & EVT_DOWN) {
    time = time*25/10;
  }
  time = step_round(time);
  set_sweep_time_us(time);
}

static void
lever_move(int status, int mode)
{
  freq_t freq = get_sweep_frequency(mode);
  if (mode == ST_SPAN){
    if (uistat.auto_center_marker) {
      freq = get_marker_frequency(active_marker);
      search_maximum(active_marker, freq, 10 );
      if (freq == 0) return;
      set_sweep_frequency(ST_CENTER, freq);
      return;
    }
    if (status & EVT_UP  ) freq = setting.frequency_var ? (freq + setting.frequency_var) : step_round(freq*4 + 1);
    if (status & EVT_DOWN) freq = setting.frequency_var ? (freq - setting.frequency_var) : step_round(freq   - 1);
  }
  else {
    freq_t span = setting.frequency_var ? setting.frequency_var : step_round(get_sweep_frequency(ST_SPAN) / 4);
    if (status & EVT_UP  ) freq+= span;
    if (status & EVT_DOWN) freq-= span;
  }
  if (freq > STOP_MAX || freq < START_MIN) return;
  set_sweep_frequency(mode, freq);
}

#define STEPRATIO 0.2
static void
ui_process_normal_lever(void)
{
  int status = btn_check();
  if (status != 0) {
    if (status & EVT_BUTTON_SINGLE_CLICK) {
      ui_mode_menu();
    } else {
      switch (uistat.lever_mode) {
      case LM_MARKER: lever_move_marker(status);   break;
      case LM_SEARCH: lever_search_marker(status); break;
      case LM_CENTER: lever_move(status, FREQ_IS_STARTSTOP() ? ST_START : ST_CENTER); break;
      case LM_SPAN:
        if (FREQ_IS_CW()) 
          lever_zoom_time(status);
        else
          lever_move(status, FREQ_IS_STARTSTOP() ? ST_STOP : ST_SPAN);
        break;
      }
    }
  }
}

#ifdef __LISTEN__
bool
ui_process_listen_lever(void)
{
  int status = btn_check();
  if (status != 0) {
    if (status & EVT_BUTTON_SINGLE_CLICK) {
      return false;
    } else {
      lever_move_marker(status);
    }
  }
  return true;
}
#endif

static void
ui_process_menu_lever(void)
{
  // Flag show, can close menu if user come out from it
  // if false user must select some thing
  const menuitem_t *menu = menu_stack[menu_current_level];
  int status = btn_check();
  if (status == 0) return;
  if (selection >=0 && status & EVT_BUTTON_SINGLE_CLICK) {
    menu_invoke(selection);
    return;
  }
  uint16_t count = current_menu_get_count();
  do {
    uint32_t mask = 1<<selection;
    if (status & EVT_UP  ) selection++;
    if (status & EVT_DOWN) selection--;
    // not close if type = form menu
    if ((uint16_t)selection >= count && !(menu[0].type & MT_FORM)){
      ui_mode_normal();
      return;
    }
    ensure_selection();
    draw_menu_mask(mask|(1<<selection));
    chThdSleepMilliseconds(100); // Add delay for not move so fast in menu
  } while ((status = btn_wait_release()) != 0);
  return;
}

static int
keypad_click(int key)
{
  int c = keypads[key].c;
  if ((c >= KP_X1 && c <= KP_G) || c == KP_m || c == KP_u || c == KP_n) {
#if 0
    float scale = 1.0;
    if (c >= KP_X1 && c <= KP_G) {
      int n = c - KP_X1;
      while (n-- > 0)
        scale *= 1000.0;
    } else if (c == KP_m) {
      scale /= 1000.0;
    } else if (c == KP_u) {
      scale /= 1000000.0;
    } else if (c == KP_n) {
      scale /= 1000000000.0;
    }
    /* numeric input done */
    uistat.value = my_atof(kp_buf) * scale;
#else
    char modifier = 0;
    if (c == KP_K) modifier = 'k';
    else if (c == KP_M) modifier = 'M';
    else if (c == KP_G) modifier = 'G';
    else if (c == KP_m) modifier = 'm';
    else if (c == KP_u) modifier = 'u';
    else if (c == KP_n) modifier = 'n';
    if (modifier) kp_buf[kp_index++] = modifier;
    kp_buf[kp_index++] = 0;
    uistat.value = my_atof(kp_buf);
    uistat.freq_value = my_atoui(kp_buf);
#endif
    set_numeric_value();
    return KP_DONE;
  } else if (c <= 9 && kp_index < NUMINPUT_LEN) {
    kp_buf[kp_index++] = '0' + c;
  } else if (c>=KP_0) {
    kp_buf[kp_index++] = keypad_scale_text[c-KP_0][0];
    if (c >=KP_10)
      kp_buf[kp_index++] = '0';
    if (c >=KP_100)
      kp_buf[kp_index++] = '0';
  } else if (c == KP_PERIOD && kp_index < NUMINPUT_LEN) {
    // check period in former input
    int j;
    for (j = 0; j < kp_index && kp_buf[j] != '.'; j++)
      ;
    // append period if there are no period
    if (kp_index == j)
      kp_buf[kp_index++] = '.';
  } else if (c == KP_MINUS && kp_index < NUMINPUT_LEN) {
    if (kp_index == 0)
      kp_buf[kp_index++] = '-';
    else {
      // always allow sign change, even when not on first position
      if (kp_buf[0] == '-') {
        kp_index = 0;
        do {
          kp_buf[kp_index] = kp_buf[kp_index+1];
          kp_index++;
        } while (kp_buf[kp_index]);
      } else {
        int j = kp_index;
        do {
          kp_buf[j+1] = kp_buf[j];
          j--;
        } while (j >= 0);
        kp_buf[0] = '-';
        kp_index++;
      }

    }
  } else if (c == KP_BS) {
    if (kp_index == 0) {
      return KP_CANCEL;
    }
    --kp_index;
  }
  kp_buf[kp_index] = '\0';
  draw_numeric_input(kp_buf);
  return KP_CONTINUE;
}

static int
keypad_apply_touch(void)
{
  int touch_x, touch_y;
  int i = 0;

  touch_position(&touch_x, &touch_y);

  while (keypads[i].c >= 0) {
    int x = KP_GET_X(keypads[i].x);
    int y = KP_GET_Y(keypads[i].y);
    if (x < touch_x && touch_x < x+KP_WIDTH && y < touch_y && touch_y < y+KP_HEIGHT) {
      uint32_t mask = (1<<i)|(1<<selection);
      // draw focus
      selection = i;
      draw_keypad(mask);
      touch_wait_release();
      // erase focus
      selection = -1;
      draw_keypad((1<<i));
      return i;
    }
    i++;
  }
  return -1;
}

static void
ui_process_keypad(void)
{
  int status;
  kp_index = 0;
  int keypads_last_index;
  for (keypads_last_index = 0; keypads[keypads_last_index+1].c >= 0; keypads_last_index++)
    ;
  while (TRUE) {
    status = btn_check();
    if (status & (EVT_UP|EVT_DOWN)) {
      int s = status;
      do {
        uint32_t mask = (1<<selection);
        if (s & EVT_UP)
          if (--selection < 0)
            selection = keypads_last_index;
        if (s & EVT_DOWN)
          if (++selection > keypads_last_index)
            selection = 0;
        draw_keypad(mask|(1<<selection));
        chThdSleepMilliseconds(100);
      } while ((s = btn_wait_release()) != 0);
    }

    if (status == EVT_BUTTON_SINGLE_CLICK) {
      if (keypad_click(selection))
        /* exit loop on done or cancel */
        break;
    }

    if (touch_check() == EVT_TOUCH_PRESSED) {
      int key = keypad_apply_touch();
      if (key >= 0 && keypad_click(key))
        /* exit loop on done or cancel */
        break;
    }
  }
  kp_help_text = NULL;
  redraw_frame();
  if (current_menu_is_form()) {
    ui_mode_menu(); //Reactivate menu after keypad
    selection = -1;
  } else {
    ui_mode_normal();
  }
  //redraw_all();
}

static void
ui_process_lever(void)
{
  switch (ui_mode) {
  case UI_NORMAL:
    ui_process_normal_lever();
    break;
  case UI_MENU:
    ui_process_menu_lever();
    break;
//  case UI_KEYPAD:
//    ui_process_keypad();
//    break;
  }
}

static void
drag_marker(int t, int m)
{
  /* wait touch release */
  do {
    int touch_x, touch_y;
    int index;
    touch_position(&touch_x, &touch_y);
    touch_x -= OFFSETX;
    touch_y -= OFFSETY;
    index = search_nearest_index(touch_x, touch_y, t);
    if (index >= 0) {
      markers[m].index = index;
      markers[m].frequency = getFrequency(index);
      redraw_marker(m);
    }
  } while (touch_check()!= EVT_TOUCH_RELEASED);
}

static int
touch_pickup_marker(int touch_x, int touch_y)
{
  int m, t;
  touch_x -= OFFSETX;
  touch_y -= OFFSETY;

  int i = MARKER_INVALID, mt;
  int min_dist = MARKER_PICKUP_DISTANCE * MARKER_PICKUP_DISTANCE;
  // Search closest marker to touch position
  for (t = 0; t < TRACES_MAX; t++) {
    if (IS_TRACE_DISABLE(t))
      continue;
    for (m = 0; m < MARKERS_MAX; m++) {
      if (!markers[m].enabled)
        continue;
      // Get distance to marker from touch point
      int dist = distance_to_index(t, markers[m].index, touch_x, touch_y);
      if (dist < min_dist) {
        min_dist = dist;
        i  = m;
        mt = t;
      }
    }
  }
  // Marker not found
  if (i == MARKER_INVALID)
    return FALSE;
  // Marker found, set as active and start drag it
  if (active_marker != i) {
    previous_marker = active_marker;
    active_marker = i;
  }
  // Disable tracking
  markers[i].mtype &= ~M_TRACKING;    // Disable tracking when dragging marker
  // Leveler mode = marker move
  select_lever_mode(LM_MARKER);
  // select trace
  uistat.current_trace = mt;
  // drag marker until release
  drag_marker(mt, i);
  return TRUE;
}

static int touch_quick_menu(int touch_x, int touch_y)
{
  if (touch_x <OFFSETX)
  {
    touch_wait_release();
    return invoke_quick_menu(touch_y);
  }
  return FALSE;
}

#ifdef __USE_SD_CARD__
//*******************************************************************************************
// Bitmap file header for LCD_WIDTH x LCD_HEIGHT image 16bpp (v4 format allow set RGB mask)
//*******************************************************************************************
#define BMP_UINT32(val)  ((val)>>0)&0xFF, ((val)>>8)&0xFF, ((val)>>16)&0xFF, ((val)>>24)&0xFF
#define BMP_H1_SIZE      (14)                        // BMP header 14 bytes
#define BMP_V4_SIZE      (56)                        // v4  header 56 bytes
#define BMP_HEAD_SIZE    (BMP_H1_SIZE + BMP_V4_SIZE) // Size of all headers
#define BMP_SIZE         (2*LCD_WIDTH*LCD_HEIGHT)    // Bitmap size = 2*w*h
#define BMP_FILE_SIZE    (BMP_SIZE + BMP_HEAD_SIZE)  // File size = headers + bitmap
static const uint8_t bmp_header_v4[14+56] = {
// BITMAPFILEHEADER (14 byte size)
  0x42, 0x4D,                // BM signature
  BMP_UINT32(BMP_FILE_SIZE), // File size (h + v4 + bitmap)
  0x00, 0x00,                // reserved
  0x00, 0x00,                // reserved
  BMP_UINT32(BMP_HEAD_SIZE), // Size of all headers (h + v4)
// BITMAPINFOv4 (56 byte size)
  BMP_UINT32(BMP_V4_SIZE),   // Data offset after this point (v4 size)
  BMP_UINT32(LCD_WIDTH),     // Width
  BMP_UINT32(LCD_HEIGHT),    // Height
  0x01, 0x00,                // Planes
  0x10, 0x00,                // 16bpp
  0x03, 0x00, 0x00, 0x00,    // Compression (BI_BITFIELDS)
  BMP_UINT32(BMP_SIZE),      // Bitmap size (w*h*2)
  0xC4, 0x0E, 0x00, 0x00,    // x Resolution (96 DPI = 96 * 39.3701 inches per metre = 0x0EC4)
  0xC4, 0x0E, 0x00, 0x00,    // y Resolution (96 DPI = 96 * 39.3701 inches per metre = 0x0EC4)
  0x00, 0x00, 0x00, 0x00,    // Palette size
  0x00, 0x00, 0x00, 0x00,    // Palette used
// Extend v4 header data (color mask for RGB565)
  0x00, 0xF8, 0x00, 0x00,    // R mask = 0b11111000 00000000
  0xE0, 0x07, 0x00, 0x00,    // G mask = 0b00000111 11100000
  0x1F, 0x00, 0x00, 0x00,    // B mask = 0b00000000 00011111
  0x00, 0x00, 0x00, 0x00     // A mask = 0b00000000 00000000
};

FRESULT open_file(char *ext)
{
#ifdef __DISABLE_HOT_INSERT__
  if (!sd_card_inserted_at_boot) {
    drawMessageBox("Warning:", "Restart tinySA to use SD card", 2000);
    return FR_NOT_READY;
  }
#endif
  FRESULT res = f_mount(fs_volume, "", 1);
  // fs_volume, fs_file and fs_filename stored at end of spi_buffer!!!!!
//  shell_printf("Mount = %d\r\n", res);
  if (res != FR_OK)
    return res;
#if FF_USE_LFN >= 1
  uint32_t tr = rtc_get_tr_bcd(); // TR read first
  uint32_t dr = rtc_get_dr_bcd(); // DR read second
  plot_printf(fs_filename, FF_LFN_BUF, "SA_%06x_%06x.%s", dr, tr, ext);
#else
  plot_printf(fs_filename, FF_LFN_BUF, "%08x.%s", rtc_get_FAT(), ext);
#endif
  res = f_open(fs_file, fs_filename, FA_CREATE_ALWAYS | FA_READ | FA_WRITE);
 return res;
}

void close_file(FRESULT res)
{
  if (res == FR_OK)
    res = f_close(fs_file);
//  time = chVTGetSystemTimeX() - time;
//  shell_printf("Total time: %dms (write %d byte/sec)\r\n", time/10, (LCD_WIDTH*LCD_HEIGHT*sizeof(uint16_t)+sizeof(bmp_header_v4))*10000/time);
  drawMessageBox("Save:", res == FR_OK ? fs_filename : "  Write failed  ", 2000);
  redraw_request|= REDRAW_AREA;
}
static bool
made_screenshot(int touch_x, int touch_y)
{
  int y, i;
  UINT size;
  if (touch_y < SD_CARD_START || touch_y > SD_CARD_START + 20 || touch_x > OFFSETX)
    return FALSE;
  ili9341_set_background(LCD_BG_COLOR);
  ili9341_fill(4, SD_CARD_START, 16, 16);
  touch_wait_release();
#ifdef __DISABLE_HOT_INSERT__
  if (!sd_card_inserted_at_boot) {
    drawMessageBox("Warning:", "Restart tinySA to use SD card", 2000);
    return FALSE;
  }
#endif
//  uint32_t time = chVTGetSystemTimeX();
//  shell_printf("Screenshot\r\n");
  FRESULT res = open_file("bmp");
  uint16_t *buf = (uint16_t *)spi_buffer;
//  shell_printf("Open %s, result = %d\r\n", fs_filename, res);
  if (res == FR_OK){
    res = f_write(fs_file, bmp_header_v4, sizeof(bmp_header_v4), &size);
    for (y = LCD_HEIGHT-1; y >= 0 && res == FR_OK; y--) {
      ili9341_read_memory(0, y, LCD_WIDTH, 1, buf);
      for (i = 0; i < LCD_WIDTH; i++)
        buf[i] = __REVSH(buf[i]); // swap byte order (example 0x10FF to 0xFF10)
      res = f_write(fs_file, buf, LCD_WIDTH*sizeof(uint16_t), &size);
    }
//   res = f_close(fs_file);
//    shell_printf("Close %d\r\n", res);
//    testLog();
  }
  close_file(res);
  return TRUE;
}

void save_to_sd(int mask)
{
  FRESULT res = open_file("csv");
  UINT size;
  if (res == FR_OK) {
    for (int i = 0; i < sweep_points; i++) {
      char *buf = (char *)spi_buffer;
      if (mask & 1) buf += plot_printf(buf, 100, "%U, ", getFrequency(i));
      if (mask & 2) buf += plot_printf(buf, 100, "%f ", value(measured[TRACE_ACTUAL][i]));
      if (mask & 4) buf += plot_printf(buf, 100, "%f ", value(measured[TRACE_STORED][i]));
      if (mask & 8) buf += plot_printf(buf, 100, "%f", value(measured[TRACE_TEMP][i]));
      buf += plot_printf(buf, 100, "\r\n");
      res = f_write(fs_file, (char *)spi_buffer, buf - (char *)spi_buffer, &size);
      if (res != FR_OK)
        break;
    }
  }
  close_file(res);
}

#endif

static int
touch_lever_mode_select(int touch_x, int touch_y)
{
  if (touch_y > HEIGHT) {
    touch_wait_release();
    // Touch on left frequency field side
    if (touch_x < FREQUENCIES_XPOS2 - 50) {
      if (uistat.lever_mode == LM_CENTER){
        ui_mode_keypad(FREQ_IS_CENTERSPAN() ? KM_CENTER : KM_START);
        return TRUE;
      }
    }
    else if (touch_x <  FREQUENCIES_XPOS2 + 50) {
      // toggle frequency mode start/stop <=> center/span
      setting.freq_mode^= FREQ_MODE_CENTER_SPAN;
      redraw_request |= REDRAW_FREQUENCY;
      return true;
    }
    else if (uistat.lever_mode == LM_SPAN) {
      ui_mode_keypad(FREQ_IS_CW() ? KM_SWEEP_TIME : (FREQ_IS_CENTERSPAN() ? KM_SPAN : KM_STOP));
      return TRUE;
    }
    select_lever_mode(touch_x < FREQUENCIES_XPOS2 ? LM_CENTER : LM_SPAN);
    return TRUE;
  }
  if (touch_x < OFFSETX)
  {
    return invoke_quick_menu(touch_y);
  }
  else if (touch_y < 25) {
#ifdef __VNA__
    if (touch_x < FREQUENCIES_XPOS2 && get_electrical_delay() != 0.0) {
      select_lever_mode(LM_EDELAY);
    } else {
#endif
      select_lever_mode(LM_MARKER);
#ifdef __VNA__
      }
#endif
    return TRUE;
  }
  return FALSE;
}

static int
touch_marker_select(int touch_x, int touch_y)
{
  int selected_marker = 0;
  if (current_menu_is_form() || touch_x > LCD_WIDTH-MENU_BUTTON_WIDTH || touch_x < 25 || touch_y > 30)
    return FALSE;
  if (touch_y > 15)
    selected_marker = 2;
  selected_marker += (touch_x >150 ? 1 : 0);
  for (int i = 0; i < MARKERS_MAX; i++) {
    if (markers[i].enabled) {
      if (selected_marker == 0) {
        if (active_marker == i) {
          extern const menuitem_t menu_marker_modify[];
          touch_wait_release();
          selection = -1;
          menu_current_level = 0;
          menu_push_submenu(menu_marker_modify);
          break;
        }
        active_marker = i;
        redraw_marker(active_marker);
        break;
      }
      selected_marker --;
    }
  }
  if (touch_y < 25) {
#ifdef __VNA__
    if (touch_x < FREQUENCIES_XPOS2 && get_electrical_delay() != 0.0) {
      select_lever_mode(LM_EDELAY);
    } else {
#endif
      select_lever_mode(LM_MARKER);
#ifdef __VNA__
      }
#endif
    return TRUE;
  }
  return FALSE;
}

static
void ui_process_touch(void)
{
  int touch_x, touch_y;
  int status = touch_check();
  if (status == EVT_TOUCH_PRESSED || status == EVT_TOUCH_DOWN) {
    touch_position(&touch_x, &touch_y);
    switch (ui_mode) {
    case UI_NORMAL:
#ifdef __USE_SD_CARD__
      if (made_screenshot(touch_x, touch_y))
        break;
#endif
      if (touch_quick_menu(touch_x, touch_y))
        break;
      // Try drag marker
      if (touch_pickup_marker(touch_x, touch_y))
        break;
      if (touch_marker_select(touch_x, touch_y))
        break;
      // Try select lever mode (top and bottom screen)
      if (touch_lever_mode_select(touch_x, touch_y)) {
//        touch_wait_release();
        break;
      }

      // switch menu mode after release
      touch_wait_release();
      selection = -1; // hide keyboard mode selection
      ui_mode_menu();
      break;
    case UI_MENU:
      menu_apply_touch(touch_x, touch_y);
      break;
    }
  }
}

static uint16_t previous_button_state = 0;

void
ui_process(void)
{
  int button_state = READ_PORT() & BUTTON_MASK;
  if (ui_mode == UI_NORMAL && current_menu_is_form()) {     //   Force into menu mode
    selection = -1; // hide keyboard mode selection
    ui_mode_menu();
  }
  if ((operation_requested&OP_LEVER) || previous_button_state != button_state) {
    ui_process_lever();
    previous_button_state = button_state;
    operation_requested = OP_NONE;
  }
  if (operation_requested&OP_TOUCH) {
   ui_process_touch();
    operation_requested = OP_NONE;
  }
  touch_start_watchdog();
}

/* Triggered when the button is pressed or released. The LED4 is set to ON.*/
static void extcb1(EXTDriver *extp, expchannel_t channel)
{
  (void)extp;
  (void)channel;
#ifdef  __USE_SD_CARD__
  if (channel == 12)
    SD_PowerOff();
#endif
  if (channel < 9)
    operation_requested|=OP_LEVER;
  // cur_button = READ_PORT() & BUTTON_MASK;
}

static const EXTConfig extcfg = {
  {
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_RISING_EDGE | EXT_CH_MODE_AUTOSTART | EXT_MODE_GPIOA, extcb1},
    {EXT_CH_MODE_RISING_EDGE | EXT_CH_MODE_AUTOSTART | EXT_MODE_GPIOA, extcb1},
    {EXT_CH_MODE_RISING_EDGE | EXT_CH_MODE_AUTOSTART | EXT_MODE_GPIOA, extcb1},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
#ifdef __WAIT_CTS_WHILE_SLEEPING__
    {EXT_CH_MODE_RISING_EDGE | EXT_CH_MODE_AUTOSTART | EXT_MODE_GPIOB, extcb1},
#else
    {EXT_CH_MODE_DISABLED, NULL},
#endif
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
#ifdef  __USE_SD_CARD__
    {EXT_CH_MODE_RISING_EDGE | EXT_CH_MODE_AUTOSTART | EXT_MODE_GPIOB, extcb1},
#else
    {EXT_CH_MODE_DISABLED, NULL},
#endif
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL},
    {EXT_CH_MODE_DISABLED, NULL}
  }
};

void
handle_touch_interrupt(void)
{
  operation_requested|= OP_TOUCH;
}

void
ui_init()
{
  adc_init();
  // Activates the EXT driver 1.
  extStart(&EXTD1, &extcfg);
  // Init touch subsystem
  touch_init();
}

void wait_user(void)
{
  touch_wait_released();
#if 0
  operation_requested = OP_NONE;
  while (true) {
    if (operation_requested & OP_TOUCH)
      break;
    if (operation_requested & OP_LEVER)
      break;
  }
#endif
}

int check_touched(void)
{
  int touched = false;
  if (touch_check() == EVT_TOUCH_RELEASED)
    touched = true;
  return touched;
}



#pragma GCC pop_options
