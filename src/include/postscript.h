#ifndef __BIOLIB__POSTSCRIPT_H__
#define __BIOLIB__POSTSCRIPT_H__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

unsigned char embed_nibble(unsigned char msnibl, unsigned char lsnibl){
      if(msnibl >=16 || lsnibl >= 16){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Nibble value more than 16. %d, %d\n", __func__, __FILE__, __LINE__, msnibl, lsnibl);
	    exit(EXIT_FAILURE);
      }
      return msnibl * 16 + lsnibl;

}
void tobinary(unsigned char* bin, unsigned char* srcint, int len){
      for(int i=0; i<len; ++i){
	    bin[i] = embed_nibble(srcint[2*i], srcint[2*i+1]);

      }
}

void toascii85(char* asc85, int* len85, const unsigned char* bin,  const int lenbin ) //#function 
{
      unsigned long map = 0;
      for(int i=0; i<lenbin; ++i){
	    map = map * 256 + bin[i];
      }
      for(int i=lenbin; i<4; ++i){ // zero padding.
	    map = map * 256 + 0;
      }
      if(map == 0){
	    asc85[0] = 'z';
	    *len85 = 1;
      }else{
	    for(int i=4; i>=0; --i){
		  int rem = map % 85;
		  asc85[i] = rem + 33;
		  map = (map - rem) / 85; 
	    }
	    *len85 = 5;
      }
}

void ps_heatmap_head(FILE* fp, char* accn, int nres, char* chain, int nmr_mdl)
{
      fprintf(fp, "%%!PS-Adobe-2.0\n");
      fprintf(fp, "%%%%Title: %s.ps\n", accn);
      fprintf(fp, "%%%%Creator: gnuplot 5.2 patchlevel 8\n");
      fprintf(fp, "%%%%CreationDate: Sun Apr  4 22:25:51 2021\n");
      fprintf(fp, "%%%%DocumentFonts: (atend)\n");
      fprintf(fp, "%%%%BoundingBox: 50 50 554 770\n");
      fprintf(fp, "%%%%Orientation: Landscape\n");
      fprintf(fp, "%%%%Pages: (atend)\n");
      fprintf(fp, "%%%%EndComments\n");
      fprintf(fp, "%%%%BeginProlog\n");
      fprintf(fp, "/gnudict 256 dict def\n");
      fprintf(fp, "gnudict begin\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%% The following true/false flags may be edited by hand if desired.\n");
      fprintf(fp, "%% The unit line width and grayscale image gamma correction may also be changed.\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/Color false def\n");
      fprintf(fp, "/Blacktext false def\n");
      fprintf(fp, "/Solid false def\n");
      fprintf(fp, "/Dashlength 1 def\n");
      fprintf(fp, "/Landscape true def\n");
      fprintf(fp, "/Level1 false def\n");
      fprintf(fp, "/Level3 false def\n");
      fprintf(fp, "/Rounded false def\n");
      fprintf(fp, "/ClipToBoundingBox false def\n");
      fprintf(fp, "/SuppressPDFMark false def\n");
      fprintf(fp, "/TransparentPatterns false def\n");
      fprintf(fp, "/gnulinewidth 5.000 def\n");
      fprintf(fp, "/userlinewidth gnulinewidth def\n");
      fprintf(fp, "/Gamma 1.0 def\n");
      fprintf(fp, "/BackgroundColor {-1.000 -1.000 -1.000} def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/vshift -46 def\n");
      fprintf(fp, "/dl1 {\n");
      fprintf(fp, "  10.0 Dashlength userlinewidth gnulinewidth div mul mul mul\n");
      fprintf(fp, "  Rounded { currentlinewidth 0.75 mul sub dup 0 le { pop 0.01 } if } if\n");
      fprintf(fp, "} def\n");
      fprintf(fp, "/dl2 {\n");
      fprintf(fp, "  10.0 Dashlength userlinewidth gnulinewidth div mul mul mul\n");
      fprintf(fp, "  Rounded { currentlinewidth 0.75 mul add } if\n");
      fprintf(fp, "} def\n");
      fprintf(fp, "/hpt_ 31.5 def\n");
      fprintf(fp, "/vpt_ 31.5 def\n");
      fprintf(fp, "/hpt hpt_ def\n");
      fprintf(fp, "/vpt vpt_ def\n");
      fprintf(fp, "/doclip {\n");
      fprintf(fp, "  ClipToBoundingBox {\n");
      fprintf(fp, "    newpath 50 50 moveto 554 50 lineto 554 770 lineto 50 770 lineto closepath\n");
      fprintf(fp, "    clip\n");
      fprintf(fp, "  } if\n");
      fprintf(fp, "} def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%% Gnuplot Prolog Version 5.2 (Dec 2017)\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%%/SuppressPDFMark true def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/M {moveto} bind def\n");
      fprintf(fp, "/L {lineto} bind def\n");
      fprintf(fp, "/R {rmoveto} bind def\n");
      fprintf(fp, "/V {rlineto} bind def\n");
      fprintf(fp, "/N {newpath moveto} bind def\n");
      fprintf(fp, "/Z {closepath} bind def\n");
      fprintf(fp, "/C {setrgbcolor} bind def\n");
      fprintf(fp, "/f {rlineto fill} bind def\n");
      fprintf(fp, "/g {setgray} bind def\n");
      fprintf(fp, "/Gshow {show} def   %% May be redefined later in the file to support UTF-8\n");
      fprintf(fp, "/vpt2 vpt 2 mul def\n");
      fprintf(fp, "/hpt2 hpt 2 mul def\n");
      fprintf(fp, "/Lshow {currentpoint stroke M 0 vshift R \n");
      fprintf(fp, "	Blacktext {gsave 0 setgray textshow grestore} {textshow} ifelse} def\n");
      fprintf(fp, "/Rshow {currentpoint stroke M dup stringwidth pop neg vshift R\n");
      fprintf(fp, "	Blacktext {gsave 0 setgray textshow grestore} {textshow} ifelse} def\n");
      fprintf(fp, "/Cshow {currentpoint stroke M dup stringwidth pop -2 div vshift R \n");
      fprintf(fp, "	Blacktext {gsave 0 setgray textshow grestore} {textshow} ifelse} def\n");
      fprintf(fp, "/UP {dup vpt_ mul /vpt exch def hpt_ mul /hpt exch def\n");
      fprintf(fp, "  /hpt2 hpt 2 mul def /vpt2 vpt 2 mul def} def\n");
      fprintf(fp, "/DL {Color {setrgbcolor Solid {pop []} if 0 setdash}\n");
      fprintf(fp, " {pop pop pop 0 setgray Solid {pop []} if 0 setdash} ifelse} def\n");
      fprintf(fp, "/BL {stroke userlinewidth 2 mul setlinewidth\n");
      fprintf(fp, "	Rounded {1 setlinejoin 1 setlinecap} if} def\n");
      fprintf(fp, "/AL {stroke userlinewidth 2 div setlinewidth\n");
      fprintf(fp, "	Rounded {1 setlinejoin 1 setlinecap} if} def\n");
      fprintf(fp, "/UL {dup gnulinewidth mul /userlinewidth exch def\n");
      fprintf(fp, "	dup 1 lt {pop 1} if 10 mul /udl exch def} def\n");
      fprintf(fp, "/PL {stroke userlinewidth setlinewidth\n");
      fprintf(fp, "	Rounded {1 setlinejoin 1 setlinecap} if} def\n");
      fprintf(fp, "3.8 setmiterlimit\n");
      fprintf(fp, "%% Classic Line colors (version 5.0)\n");
      fprintf(fp, "/LCw {1 1 1} def\n");
      fprintf(fp, "/LCb {0 0 0} def\n");
      fprintf(fp, "/LCa {0 0 0} def\n");
      fprintf(fp, "/LC0 {1 0 0} def\n");
      fprintf(fp, "/LC1 {0 1 0} def\n");
      fprintf(fp, "/LC2 {0 0 1} def\n");
      fprintf(fp, "/LC3 {1 0 1} def\n");
      fprintf(fp, "/LC4 {0 1 1} def\n");
      fprintf(fp, "/LC5 {1 1 0} def\n");
      fprintf(fp, "/LC6 {0 0 0} def\n");
      fprintf(fp, "/LC7 {1 0.3 0} def\n");
      fprintf(fp, "/LC8 {0.5 0.5 0.5} def\n");
      fprintf(fp, "%% Default dash patterns (version 5.0)\n");
      fprintf(fp, "/LTB {BL [] LCb DL} def\n");
      fprintf(fp, "/LTw {PL [] 1 setgray} def\n");
      fprintf(fp, "/LTb {PL [] LCb DL} def\n");
      fprintf(fp, "/LTa {AL [1 udl mul 2 udl mul] 0 setdash LCa setrgbcolor} def\n");
      fprintf(fp, "/LT0 {PL [] LC0 DL} def\n");
      fprintf(fp, "/LT1 {PL [2 dl1 3 dl2] LC1 DL} def\n");
      fprintf(fp, "/LT2 {PL [1 dl1 1.5 dl2] LC2 DL} def\n");
      fprintf(fp, "/LT3 {PL [6 dl1 2 dl2 1 dl1 2 dl2] LC3 DL} def\n");
      fprintf(fp, "/LT4 {PL [1 dl1 2 dl2 6 dl1 2 dl2 1 dl1 2 dl2] LC4 DL} def\n");
      fprintf(fp, "/LT5 {PL [4 dl1 2 dl2] LC5 DL} def\n");
      fprintf(fp, "/LT6 {PL [1.5 dl1 1.5 dl2 1.5 dl1 1.5 dl2 1.5 dl1 6 dl2] LC6 DL} def\n");
      fprintf(fp, "/LT7 {PL [3 dl1 3 dl2 1 dl1 3 dl2] LC7 DL} def\n");
      fprintf(fp, "/LT8 {PL [2 dl1 2 dl2 2 dl1 6 dl2] LC8 DL} def\n");
      fprintf(fp, "/SL {[] 0 setdash} def\n");
      fprintf(fp, "/Pnt {stroke [] 0 setdash gsave 1 setlinecap M 0 0 V stroke grestore} def\n");
      fprintf(fp, "/Dia {stroke [] 0 setdash 2 copy vpt add M\n");
      fprintf(fp, "  hpt neg vpt neg V hpt vpt neg V\n");
      fprintf(fp, "  hpt vpt V hpt neg vpt V closepath stroke\n");
      fprintf(fp, "  Pnt} def\n");
      fprintf(fp, "/Pls {stroke [] 0 setdash vpt sub M 0 vpt2 V\n");
      fprintf(fp, "  currentpoint stroke M\n");
      fprintf(fp, "  hpt neg vpt neg R hpt2 0 V stroke\n");
      fprintf(fp, " } def\n");
      fprintf(fp, "/Box {stroke [] 0 setdash 2 copy exch hpt sub exch vpt add M\n");
      fprintf(fp, "  0 vpt2 neg V hpt2 0 V 0 vpt2 V\n");
      fprintf(fp, "  hpt2 neg 0 V closepath stroke\n");
      fprintf(fp, "  Pnt} def\n");
      fprintf(fp, "/Crs {stroke [] 0 setdash exch hpt sub exch vpt add M\n");
      fprintf(fp, "  hpt2 vpt2 neg V currentpoint stroke M\n");
      fprintf(fp, "  hpt2 neg 0 R hpt2 vpt2 V stroke} def\n");
      fprintf(fp, "/TriU {stroke [] 0 setdash 2 copy vpt 1.12 mul add M\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V closepath stroke\n");
      fprintf(fp, "  Pnt} def\n");
      fprintf(fp, "/Star {2 copy Pls Crs} def\n");
      fprintf(fp, "/BoxF {stroke [] 0 setdash exch hpt sub exch vpt add M\n");
      fprintf(fp, "  0 vpt2 neg V hpt2 0 V 0 vpt2 V\n");
      fprintf(fp, "  hpt2 neg 0 V closepath fill} def\n");
      fprintf(fp, "/TriUF {stroke [] 0 setdash vpt 1.12 mul add M\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V closepath fill} def\n");
      fprintf(fp, "/TriD {stroke [] 0 setdash 2 copy vpt 1.12 mul sub M\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V closepath stroke\n");
      fprintf(fp, "  Pnt} def\n");
      fprintf(fp, "/TriDF {stroke [] 0 setdash vpt 1.12 mul sub M\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V closepath fill} def\n");
      fprintf(fp, "/DiaF {stroke [] 0 setdash vpt add M\n");
      fprintf(fp, "  hpt neg vpt neg V hpt vpt neg V\n");
      fprintf(fp, "  hpt vpt V hpt neg vpt V closepath fill} def\n");
      fprintf(fp, "/Pent {stroke [] 0 setdash 2 copy gsave\n");
      fprintf(fp, "  translate 0 hpt M 4 {72 rotate 0 hpt L} repeat\n");
      fprintf(fp, "  closepath stroke grestore Pnt} def\n");
      fprintf(fp, "/PentF {stroke [] 0 setdash gsave\n");
      fprintf(fp, "  translate 0 hpt M 4 {72 rotate 0 hpt L} repeat\n");
      fprintf(fp, "  closepath fill grestore} def\n");
      fprintf(fp, "/Circle {stroke [] 0 setdash 2 copy\n");
      fprintf(fp, "  hpt 0 360 arc stroke Pnt} def\n");
      fprintf(fp, "/CircleF {stroke [] 0 setdash hpt 0 360 arc fill} def\n");
      fprintf(fp, "/C0 {BL [] 0 setdash 2 copy moveto vpt 90 450 arc} bind def\n");
      fprintf(fp, "/C1 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 0 90 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C2 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 90 180 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C3 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 0 180 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C4 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 180 270 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C5 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 0 90 arc\n");
      fprintf(fp, "	2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 180 270 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc} bind def\n");
      fprintf(fp, "/C6 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 90 270 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C7 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 0 270 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C8 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 270 360 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C9 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 270 450 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C10 {BL [] 0 setdash 2 copy 2 copy moveto vpt 270 360 arc closepath fill\n");
      fprintf(fp, "	2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 90 180 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C11 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 0 180 arc closepath fill\n");
      fprintf(fp, "	2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 270 360 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C12 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 180 360 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C13 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 0 90 arc closepath fill\n");
      fprintf(fp, "	2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 180 360 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/C14 {BL [] 0 setdash 2 copy moveto\n");
      fprintf(fp, "	2 copy vpt 90 360 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc} bind def\n");
      fprintf(fp, "/C15 {BL [] 0 setdash 2 copy vpt 0 360 arc closepath fill\n");
      fprintf(fp, "	vpt 0 360 arc closepath} bind def\n");
      fprintf(fp, "/Rec {newpath 4 2 roll moveto 1 index 0 rlineto 0 exch rlineto\n");
      fprintf(fp, "	neg 0 rlineto closepath} bind def\n");
      fprintf(fp, "/Square {dup Rec} bind def\n");
      fprintf(fp, "/Bsquare {vpt sub exch vpt sub exch vpt2 Square} bind def\n");
      fprintf(fp, "/S0 {BL [] 0 setdash 2 copy moveto 0 vpt rlineto BL Bsquare} bind def\n");
      fprintf(fp, "/S1 {BL [] 0 setdash 2 copy vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S2 {BL [] 0 setdash 2 copy exch vpt sub exch vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S3 {BL [] 0 setdash 2 copy exch vpt sub exch vpt2 vpt Rec fill Bsquare} bind def\n");
      fprintf(fp, "/S4 {BL [] 0 setdash 2 copy exch vpt sub exch vpt sub vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S5 {BL [] 0 setdash 2 copy 2 copy vpt Square fill\n");
      fprintf(fp, "	exch vpt sub exch vpt sub vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S6 {BL [] 0 setdash 2 copy exch vpt sub exch vpt sub vpt vpt2 Rec fill Bsquare} bind def\n");
      fprintf(fp, "/S7 {BL [] 0 setdash 2 copy exch vpt sub exch vpt sub vpt vpt2 Rec fill\n");
      fprintf(fp, "	2 copy vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S8 {BL [] 0 setdash 2 copy vpt sub vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S9 {BL [] 0 setdash 2 copy vpt sub vpt vpt2 Rec fill Bsquare} bind def\n");
      fprintf(fp, "/S10 {BL [] 0 setdash 2 copy vpt sub vpt Square fill 2 copy exch vpt sub exch vpt Square fill\n");
      fprintf(fp, "	Bsquare} bind def\n");
      fprintf(fp, "/S11 {BL [] 0 setdash 2 copy vpt sub vpt Square fill 2 copy exch vpt sub exch vpt2 vpt Rec fill\n");
      fprintf(fp, "	Bsquare} bind def\n");
      fprintf(fp, "/S12 {BL [] 0 setdash 2 copy exch vpt sub exch vpt sub vpt2 vpt Rec fill Bsquare} bind def\n");
      fprintf(fp, "/S13 {BL [] 0 setdash 2 copy exch vpt sub exch vpt sub vpt2 vpt Rec fill\n");
      fprintf(fp, "	2 copy vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S14 {BL [] 0 setdash 2 copy exch vpt sub exch vpt sub vpt2 vpt Rec fill\n");
      fprintf(fp, "	2 copy exch vpt sub exch vpt Square fill Bsquare} bind def\n");
      fprintf(fp, "/S15 {BL [] 0 setdash 2 copy Bsquare fill Bsquare} bind def\n");
      fprintf(fp, "/D0 {gsave translate 45 rotate 0 0 S0 stroke grestore} bind def\n");
      fprintf(fp, "/D1 {gsave translate 45 rotate 0 0 S1 stroke grestore} bind def\n");
      fprintf(fp, "/D2 {gsave translate 45 rotate 0 0 S2 stroke grestore} bind def\n");
      fprintf(fp, "/D3 {gsave translate 45 rotate 0 0 S3 stroke grestore} bind def\n");
      fprintf(fp, "/D4 {gsave translate 45 rotate 0 0 S4 stroke grestore} bind def\n");
      fprintf(fp, "/D5 {gsave translate 45 rotate 0 0 S5 stroke grestore} bind def\n");
      fprintf(fp, "/D6 {gsave translate 45 rotate 0 0 S6 stroke grestore} bind def\n");
      fprintf(fp, "/D7 {gsave translate 45 rotate 0 0 S7 stroke grestore} bind def\n");
      fprintf(fp, "/D8 {gsave translate 45 rotate 0 0 S8 stroke grestore} bind def\n");
      fprintf(fp, "/D9 {gsave translate 45 rotate 0 0 S9 stroke grestore} bind def\n");
      fprintf(fp, "/D10 {gsave translate 45 rotate 0 0 S10 stroke grestore} bind def\n");
      fprintf(fp, "/D11 {gsave translate 45 rotate 0 0 S11 stroke grestore} bind def\n");
      fprintf(fp, "/D12 {gsave translate 45 rotate 0 0 S12 stroke grestore} bind def\n");
      fprintf(fp, "/D13 {gsave translate 45 rotate 0 0 S13 stroke grestore} bind def\n");
      fprintf(fp, "/D14 {gsave translate 45 rotate 0 0 S14 stroke grestore} bind def\n");
      fprintf(fp, "/D15 {gsave translate 45 rotate 0 0 S15 stroke grestore} bind def\n");
      fprintf(fp, "/DiaE {stroke [] 0 setdash vpt add M\n");
      fprintf(fp, "  hpt neg vpt neg V hpt vpt neg V\n");
      fprintf(fp, "  hpt vpt V hpt neg vpt V closepath stroke} def\n");
      fprintf(fp, "/BoxE {stroke [] 0 setdash exch hpt sub exch vpt add M\n");
      fprintf(fp, "  0 vpt2 neg V hpt2 0 V 0 vpt2 V\n");
      fprintf(fp, "  hpt2 neg 0 V closepath stroke} def\n");
      fprintf(fp, "/TriUE {stroke [] 0 setdash vpt 1.12 mul add M\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V closepath stroke} def\n");
      fprintf(fp, "/TriDE {stroke [] 0 setdash vpt 1.12 mul sub M\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V closepath stroke} def\n");
      fprintf(fp, "/PentE {stroke [] 0 setdash gsave\n");
      fprintf(fp, "  translate 0 hpt M 4 {72 rotate 0 hpt L} repeat\n");
      fprintf(fp, "  closepath stroke grestore} def\n");
      fprintf(fp, "/CircE {stroke [] 0 setdash \n");
      fprintf(fp, "  hpt 0 360 arc stroke} def\n");
      fprintf(fp, "/Opaque {gsave closepath 1 setgray fill grestore 0 setgray closepath} def\n");
      fprintf(fp, "/DiaW {stroke [] 0 setdash vpt add M\n");
      fprintf(fp, "  hpt neg vpt neg V hpt vpt neg V\n");
      fprintf(fp, "  hpt vpt V hpt neg vpt V Opaque stroke} def\n");
      fprintf(fp, "/BoxW {stroke [] 0 setdash exch hpt sub exch vpt add M\n");
      fprintf(fp, "  0 vpt2 neg V hpt2 0 V 0 vpt2 V\n");
      fprintf(fp, "  hpt2 neg 0 V Opaque stroke} def\n");
      fprintf(fp, "/TriUW {stroke [] 0 setdash vpt 1.12 mul add M\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V Opaque stroke} def\n");
      fprintf(fp, "/TriDW {stroke [] 0 setdash vpt 1.12 mul sub M\n");
      fprintf(fp, "  hpt neg vpt 1.62 mul V\n");
      fprintf(fp, "  hpt 2 mul 0 V\n");
      fprintf(fp, "  hpt neg vpt -1.62 mul V Opaque stroke} def\n");
      fprintf(fp, "/PentW {stroke [] 0 setdash gsave\n");
      fprintf(fp, "  translate 0 hpt M 4 {72 rotate 0 hpt L} repeat\n");
      fprintf(fp, "  Opaque stroke grestore} def\n");
      fprintf(fp, "/CircW {stroke [] 0 setdash \n");
      fprintf(fp, "  hpt 0 360 arc Opaque stroke} def\n");
      fprintf(fp, "/BoxFill {gsave Rec 1 setgray fill grestore} def\n");
      fprintf(fp, "/Density {\n");
      fprintf(fp, "  /Fillden exch def\n");
      fprintf(fp, "  currentrgbcolor\n");
      fprintf(fp, "  /ColB exch def /ColG exch def /ColR exch def\n");
      fprintf(fp, "  /ColR ColR Fillden mul Fillden sub 1 add def\n");
      fprintf(fp, "  /ColG ColG Fillden mul Fillden sub 1 add def\n");
      fprintf(fp, "  /ColB ColB Fillden mul Fillden sub 1 add def\n");
      fprintf(fp, "  ColR ColG ColB setrgbcolor} def\n");
      fprintf(fp, "/BoxColFill {gsave Rec PolyFill} def\n");
      fprintf(fp, "/PolyFill {gsave Density fill grestore grestore} def\n");
      fprintf(fp, "/h {rlineto rlineto rlineto closepath gsave fill grestore stroke} bind def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%% PostScript Level 1 Pattern Fill routine for rectangles\n");
      fprintf(fp, "%% Usage: x y w h s a XX PatternFill\n");
      fprintf(fp, "%%	x,y = lower left corner of box to be filled\n");
      fprintf(fp, "%%	w,h = width and height of box\n");
      fprintf(fp, "%%	  a = angle in degrees between lines and x-axis\n");
      fprintf(fp, "%%	 XX = 0/1 for no/yes cross-hatch\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/PatternFill {gsave /PFa [ 9 2 roll ] def\n");
      fprintf(fp, "  PFa 0 get PFa 2 get 2 div add PFa 1 get PFa 3 get 2 div add translate\n");
      fprintf(fp, "  PFa 2 get -2 div PFa 3 get -2 div PFa 2 get PFa 3 get Rec\n");
      fprintf(fp, "  TransparentPatterns {} {gsave 1 setgray fill grestore} ifelse\n");
      fprintf(fp, "  clip\n");
      fprintf(fp, "  currentlinewidth 0.5 mul setlinewidth\n");
      fprintf(fp, "  /PFs PFa 2 get dup mul PFa 3 get dup mul add sqrt def\n");
      fprintf(fp, "  0 0 M PFa 5 get rotate PFs -2 div dup translate\n");
      fprintf(fp, "  0 1 PFs PFa 4 get div 1 add floor cvi\n");
      fprintf(fp, "	{PFa 4 get mul 0 M 0 PFs V} for\n");
      fprintf(fp, "  0 PFa 6 get ne {\n");
      fprintf(fp, "	0 1 PFs PFa 4 get div 1 add floor cvi\n");
      fprintf(fp, "	{PFa 4 get mul 0 2 1 roll M PFs 0 V} for\n");
      fprintf(fp, " } if\n");
      fprintf(fp, "  stroke grestore} def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/languagelevel where\n");
      fprintf(fp, " {pop languagelevel} {1} ifelse\n");
      fprintf(fp, "dup 2 lt\n");
      fprintf(fp, "	{/InterpretLevel1 true def\n");
      fprintf(fp, "	 /InterpretLevel3 false def}\n");
      fprintf(fp, "	{/InterpretLevel1 Level1 def\n");
      fprintf(fp, "	 2 gt\n");
      fprintf(fp, "	    {/InterpretLevel3 Level3 def}\n");
      fprintf(fp, "	    {/InterpretLevel3 false def}\n");
      fprintf(fp, "	 ifelse }\n");
      fprintf(fp, " ifelse\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%% PostScript level 2 pattern fill definitions\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/Level2PatternFill {\n");
      fprintf(fp, "/Tile8x8 {/PaintType 2 /PatternType 1 /TilingType 1 /BBox [0 0 8 8] /XStep 8 /YStep 8}\n");
      fprintf(fp, "	bind def\n");
      fprintf(fp, "/KeepColor {currentrgbcolor [/Pattern /DeviceRGB] setcolorspace} bind def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop 0 0 M 8 8 L 0 8 M 8 0 L stroke} \n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat1 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop 0 0 M 8 8 L 0 8 M 8 0 L stroke\n");
      fprintf(fp, "	0 4 M 4 8 L 8 4 L 4 0 L 0 4 L stroke}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat2 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop 0 0 M 0 8 L\n");
      fprintf(fp, "	8 8 L 8 0 L 0 0 L fill}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat3 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop -4 8 M 8 -4 L\n");
      fprintf(fp, "	0 12 M 12 0 L stroke}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat4 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop -4 0 M 8 12 L\n");
      fprintf(fp, "	0 -4 M 12 8 L stroke}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat5 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop -2 8 M 4 -4 L\n");
      fprintf(fp, "	0 12 M 8 -4 L 4 12 M 10 0 L stroke}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat6 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop -2 0 M 4 12 L\n");
      fprintf(fp, "	0 -4 M 8 12 L 4 -4 M 10 8 L stroke}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat7 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop 8 -2 M -4 4 L\n");
      fprintf(fp, "	12 0 M -4 8 L 12 4 M 0 10 L stroke}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat8 exch def\n");
      fprintf(fp, "<< Tile8x8\n");
      fprintf(fp, " /PaintProc {0.5 setlinewidth pop 0 -2 M 12 4 L\n");
      fprintf(fp, "	-4 0 M 12 8 L -4 4 M 8 10 L stroke}\n");
      fprintf(fp, ">> matrix makepattern\n");
      fprintf(fp, "/Pat9 exch def\n");
      fprintf(fp, "/Pattern1 {PatternBgnd KeepColor Pat1 setpattern} bind def\n");
      fprintf(fp, "/Pattern2 {PatternBgnd KeepColor Pat2 setpattern} bind def\n");
      fprintf(fp, "/Pattern3 {PatternBgnd KeepColor Pat3 setpattern} bind def\n");
      fprintf(fp, "/Pattern4 {PatternBgnd KeepColor Landscape {Pat5} {Pat4} ifelse setpattern} bind def\n");
      fprintf(fp, "/Pattern5 {PatternBgnd KeepColor Landscape {Pat4} {Pat5} ifelse setpattern} bind def\n");
      fprintf(fp, "/Pattern6 {PatternBgnd KeepColor Landscape {Pat9} {Pat6} ifelse setpattern} bind def\n");
      fprintf(fp, "/Pattern7 {PatternBgnd KeepColor Landscape {Pat8} {Pat7} ifelse setpattern} bind def\n");
      fprintf(fp, "} def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%%End of PostScript Level 2 code\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/PatternBgnd {\n");
      fprintf(fp, "  TransparentPatterns {} {gsave 1 setgray fill grestore} ifelse\n");
      fprintf(fp, "} def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%% Substitute for Level 2 pattern fill codes with\n");
      fprintf(fp, "%% grayscale if Level 2 support is not selected.\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/Level1PatternFill {\n");
      fprintf(fp, "/Pattern1 {0.250 Density} bind def\n");
      fprintf(fp, "/Pattern2 {0.500 Density} bind def\n");
      fprintf(fp, "/Pattern3 {0.750 Density} bind def\n");
      fprintf(fp, "/Pattern4 {0.125 Density} bind def\n");
      fprintf(fp, "/Pattern5 {0.375 Density} bind def\n");
      fprintf(fp, "/Pattern6 {0.625 Density} bind def\n");
      fprintf(fp, "/Pattern7 {0.875 Density} bind def\n");
      fprintf(fp, "} def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%% Now test for support of Level 2 code\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "Level1 {Level1PatternFill} {Level2PatternFill} ifelse\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/Symbol-Oblique /Symbol findfont [1 0 .167 1 0 0] makefont\n");
      fprintf(fp, "dup length dict begin {1 index /FID eq {pop pop} {def} ifelse} forall\n");
      fprintf(fp, "currentdict end definefont pop\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/Metrics {ExtendTextBox Gswidth} def\n");
      fprintf(fp, "/Lwidth {currentpoint stroke M 0 vshift R Metrics} def\n");
      fprintf(fp, "/Rwidth {currentpoint stroke M dup stringwidth pop neg vshift R Metrics} def\n");
      fprintf(fp, "/Cwidth {currentpoint stroke M dup stringwidth pop -2 div vshift R Metrics} def\n");
      fprintf(fp, "/GLwidth {currentpoint stroke M 0 vshift R {ExtendTextBox} forall} def\n");
      fprintf(fp, "/GRwidth {currentpoint stroke M dup Gwidth vshift R {ExtendTextBox} forall} def\n");
      fprintf(fp, "/GCwidth {currentpoint stroke M dup Gwidth 2 div vshift R {ExtendTextBox} forall} def\n");
      fprintf(fp, "/GLwidth2 {0 Gwidth AddGlyphWidth} def\n");
      fprintf(fp, "/GRwidth2 {Gwidth -1 mul 0 AddGlyphWidth} def\n");
      fprintf(fp, "/GCwidth2 {Gwidth 2 div dup -1 mul AddGlyphWidth} def\n");
      fprintf(fp, "/AddGlyphWidth { dup TBx2 gt {userdict /TBx2 3 -1 roll put} {pop} ifelse\n");
      fprintf(fp, "                 dup TBx1 lt {userdict /TBx1 3 -1 roll put} {pop} ifelse } def\n");
      fprintf(fp, "/MFshow {\n");
      fprintf(fp, "   { dup 5 get 3 ge\n");
      fprintf(fp, "     { 5 get 3 eq {gsave} {grestore} ifelse }\n");
      fprintf(fp, "     {dup dup 0 get findfont exch 1 get scalefont setfont\n");
      fprintf(fp, "     [ currentpoint ] exch dup 2 get 0 exch R dup 5 get 2 ne {dup dup 6\n");
      fprintf(fp, "     get exch 4 get {textshow} {Metrics pop 0 R} ifelse }if dup 5 get 0 eq\n");
      fprintf(fp, "     {dup 3 get {2 get neg 0 exch R pop} {pop aload pop M} ifelse} {dup 5\n");
      fprintf(fp, "     get 1 eq {dup 2 get exch dup 3 get exch 6 get Gswidth pop -2 div\n");
      fprintf(fp, "     dup 0 R} {dup 6 get Gswidth pop -2 div 0 R 6 get\n");
      fprintf(fp, "     textshow 2 index {aload pop M neg 3 -1 roll neg R pop pop} {pop pop pop\n");
      fprintf(fp, "     pop aload pop M} ifelse }ifelse }ifelse }\n");
      fprintf(fp, "     ifelse }\n");
      fprintf(fp, "   forall} def\n");
      fprintf(fp, "/Gswidth {dup type /stringtype eq {stringwidth} {pop (n) stringwidth} ifelse} def\n");
      fprintf(fp, "/MFwidth {0 exch { dup 5 get 3 ge { 5 get 3 eq { 0 } { pop } ifelse }\n");
      fprintf(fp, " {dup 3 get{dup dup 0 get findfont exch 1 get scalefont setfont\n");
      fprintf(fp, "     6 get Gswidth pop add} {pop} ifelse} ifelse} forall} def\n");
      fprintf(fp, "/MLshow { currentpoint stroke M\n");
      fprintf(fp, "  0 exch R\n");
      fprintf(fp, "  Blacktext {gsave 0 setgray MFshow grestore} {MFshow} ifelse } bind def\n");
      fprintf(fp, "/MRshow { currentpoint stroke M\n");
      fprintf(fp, "  exch dup MFwidth neg 3 -1 roll R\n");
      fprintf(fp, "  Blacktext {gsave 0 setgray MFshow grestore} {MFshow} ifelse } bind def\n");
      fprintf(fp, "/MCshow { currentpoint stroke M\n");
      fprintf(fp, "  exch dup MFwidth -2 div 3 -1 roll R\n");
      fprintf(fp, "  Blacktext {gsave 0 setgray MFshow grestore} {MFshow} ifelse } bind def\n");
      fprintf(fp, "/XYsave    { [( ) 1 2 true false 3 ()] } bind def\n");
      fprintf(fp, "/XYrestore { [( ) 1 2 true false 4 ()] } bind def\n");
      fprintf(fp, "Level1 SuppressPDFMark or \n");
      fprintf(fp, "{} {\n");
      fprintf(fp, "/SDict 10 dict def\n");
      fprintf(fp, "systemdict /pdfmark known not {\n");
      fprintf(fp, "  userdict /pdfmark systemdict /cleartomark get put\n");
      fprintf(fp, "} if\n");
      fprintf(fp, "SDict begin [\n");
      fprintf(fp, "  /Title (%s.ps)\n", accn);
      fprintf(fp, "  /Subject (gnuplot plot)\n");
      fprintf(fp, "  /Creator (gnuplot 5.2 patchlevel 8)\n");
      fprintf(fp, "%%  /Producer (gnuplot)\n");
      fprintf(fp, "%%  /Keywords ()\n");
      fprintf(fp, "  /CreationDate (Sun Apr  4 22:25:51 2021)\n");
      fprintf(fp, "  /DOCINFO pdfmark\n");
      fprintf(fp, "end\n");
      fprintf(fp, "} ifelse\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "%% Support for boxed text - Ethan A Merritt Sep 2016\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "/InitTextBox { userdict /TBy2 3 -1 roll put userdict /TBx2 3 -1 roll put\n");
      fprintf(fp, "           userdict /TBy1 3 -1 roll put userdict /TBx1 3 -1 roll put\n");
      fprintf(fp, "	   /Boxing true def } def\n");
      fprintf(fp, "/ExtendTextBox { dup type /stringtype eq\n");
      fprintf(fp, "    { Boxing { gsave dup false charpath pathbbox\n");
      fprintf(fp, "      dup TBy2 gt {userdict /TBy2 3 -1 roll put} {pop} ifelse\n");
      fprintf(fp, "      dup TBx2 gt {userdict /TBx2 3 -1 roll put} {pop} ifelse\n");
      fprintf(fp, "      dup TBy1 lt {userdict /TBy1 3 -1 roll put} {pop} ifelse\n");
      fprintf(fp, "      dup TBx1 lt {userdict /TBx1 3 -1 roll put} {pop} ifelse\n");
      fprintf(fp, "      grestore } if }\n");
      fprintf(fp, "    {} ifelse} def\n");
      fprintf(fp, "/PopTextBox { newpath TBx1 TBxmargin sub TBy1 TBymargin sub M\n");
      fprintf(fp, "               TBx1 TBxmargin sub TBy2 TBymargin add L\n");
      fprintf(fp, "	       TBx2 TBxmargin add TBy2 TBymargin add L\n");
      fprintf(fp, "	       TBx2 TBxmargin add TBy1 TBymargin sub L closepath } def\n");
      fprintf(fp, "/DrawTextBox { PopTextBox stroke /Boxing false def} def\n");
      fprintf(fp, "/FillTextBox { gsave PopTextBox fill grestore /Boxing false def} def\n");
      fprintf(fp, "0 0 0 0 InitTextBox\n");
      fprintf(fp, "/TBxmargin 20 def\n");
      fprintf(fp, "/TBymargin 20 def\n");
      fprintf(fp, "/Boxing false def\n");
      fprintf(fp, "/textshow { ExtendTextBox Gshow } def\n");
      fprintf(fp, "%%\n");
      fprintf(fp, "end\n");
      fprintf(fp, "%%%%EndProlog\n");
      fprintf(fp, "%%%%Page: 1 1\n");
      fprintf(fp, "gnudict begin\n");
      fprintf(fp, "gsave\n");
      fprintf(fp, "doclip\n");
      fprintf(fp, "50 50 translate\n");
      fprintf(fp, "0.100 0.100 scale\n");
      fprintf(fp, "90 rotate\n");
      fprintf(fp, "0 -5040 translate\n");
      fprintf(fp, "0 setgray\n");
      fprintf(fp, "newpath\n");
      fprintf(fp, "(Helvetica) findfont 140 scalefont setfont\n");
      fprintf(fp, "BackgroundColor 0 lt 3 1 roll 0 lt exch 0 lt or or not {gsave BackgroundColor C clippath fill grestore} if\n");
      fprintf(fp, "gsave %% colour palette begin\n");
      fprintf(fp, "/maxcolors 9 def\n");
      fprintf(fp, "/HSV2RGB {  exch dup 0.0 eq {pop exch pop dup dup} %% achromatic gray\n");
      fprintf(fp, "  { /HSVs exch def /HSVv exch def 6.0 mul dup floor dup 3 1 roll sub\n");
      fprintf(fp, "     /HSVf exch def /HSVi exch cvi def /HSVp HSVv 1.0 HSVs sub mul def\n");
      fprintf(fp, "	 /HSVq HSVv 1.0 HSVs HSVf mul sub mul def \n");
      fprintf(fp, "	 /HSVt HSVv 1.0 HSVs 1.0 HSVf sub mul sub mul def\n");
      fprintf(fp, "	 /HSVi HSVi 6 mod def 0 HSVi eq {HSVv HSVt HSVp}\n");
      fprintf(fp, "	 {1 HSVi eq {HSVq HSVv HSVp}{2 HSVi eq {HSVp HSVv HSVt}\n");
      fprintf(fp, "	 {3 HSVi eq {HSVp HSVq HSVv}{4 HSVi eq {HSVt HSVp HSVv}\n");
      fprintf(fp, "	 {HSVv HSVp HSVq} ifelse} ifelse} ifelse} ifelse} ifelse\n");
      fprintf(fp, "  } ifelse} def\n");
      fprintf(fp, "/Constrain {\n");
      fprintf(fp, "  dup 0 lt {0 exch pop}{dup 1 gt {1 exch pop} if} ifelse} def\n");
      fprintf(fp, "/CMY2RGB {  1 exch sub exch 1 exch sub 3 2 roll 1 exch sub 3 1 roll exch } def\n");
      fprintf(fp, "/XYZ2RGB {  3 copy -0.9017 mul exch -0.1187 mul add exch 0.0585 mul exch add\n");
      fprintf(fp, "  Constrain 4 1 roll 3 copy -0.0279 mul exch 1.999 mul add exch\n");
      fprintf(fp, "  -0.9844 mul add Constrain 5 1 roll -0.2891 mul exch -0.5338 mul add\n");
      fprintf(fp, "  exch 1.91 mul exch add Constrain 3 1 roll} def\n");
      fprintf(fp, "/SelectSpace {ColorSpace (HSV) eq {HSV2RGB}{ColorSpace (XYZ) eq {\n");
      fprintf(fp, "  XYZ2RGB}{ColorSpace (CMY) eq {CMY2RGB}{ColorSpace (YIQ) eq {YIQ2RGB}\n");
      fprintf(fp, "  if} ifelse} ifelse} ifelse} def\n");
      fprintf(fp, "/InterpolatedColor true def\n");
      fprintf(fp, "/grayindex {/gidx 0 def\n");
      fprintf(fp, "  {GrayA gidx get grayv ge {exit} if /gidx gidx 1 add def} loop} def\n");
      fprintf(fp, "/dgdx {grayv GrayA gidx get sub GrayA gidx 1 sub get\n");
      fprintf(fp, "  GrayA gidx get sub div} def \n");
      fprintf(fp, "/redvalue {RedA gidx get RedA gidx 1 sub get\n");
      fprintf(fp, "  RedA gidx get sub dgdxval mul add} def\n");
      fprintf(fp, "/greenvalue {GreenA gidx get GreenA gidx 1 sub get\n");
      fprintf(fp, "  GreenA gidx get sub dgdxval mul add} def\n");
      fprintf(fp, "/bluevalue {BlueA gidx get BlueA gidx 1 sub get\n");
      fprintf(fp, "  BlueA gidx get sub dgdxval mul add} def\n");
      fprintf(fp, "/interpolate {\n");
      fprintf(fp, "  grayindex grayv GrayA gidx get sub abs 1e-5 le\n");
      fprintf(fp, "    {RedA gidx get GreenA gidx get BlueA gidx get}\n");
      fprintf(fp, "    {/dgdxval dgdx def redvalue greenvalue bluevalue} ifelse} def\n");
      fprintf(fp, "/GrayA [0 .125 .25 .375 .5 .625 .75 .875 1 ] def\n");
      fprintf(fp, "/RedA [.1451 1 .6588 1 1 0 0 1 0 ] def\n");
      fprintf(fp, "/GreenA [.1294 .1686 .6588 .7373 0 .3333 1 1 0 ] def\n");
      fprintf(fp, "/BlueA [.1294 .8353 .6588 0 0 1 0 1 0 ] def\n");
      fprintf(fp, "/pm3dround {maxcolors 0 gt {dup 1 ge\n");
      fprintf(fp, "	{pop 1} {maxcolors mul floor maxcolors 1 sub div} ifelse} if} def\n");
      fprintf(fp, "/pm3dGamma 1.0 1.5 Gamma mul div def\n");
      fprintf(fp, "/ColorSpace (RGB) def\n");
      fprintf(fp, "Color InterpolatedColor or { %% COLOUR vs. GRAY map\n");
      fprintf(fp, "  InterpolatedColor { %%%% Interpolation vs. RGB-Formula\n");
      fprintf(fp, "    /g {stroke pm3dround /grayv exch def interpolate\n");
      fprintf(fp, "        SelectSpace setrgbcolor} bind def\n");
      fprintf(fp, "  }{\n");
      fprintf(fp, "  /g {stroke pm3dround dup cF7 Constrain exch dup cF5 Constrain exch cF15 Constrain \n");
      fprintf(fp, "       SelectSpace setrgbcolor} bind def\n");
      fprintf(fp, "  } ifelse\n");
      fprintf(fp, "}{\n");
      fprintf(fp, "  /g {stroke pm3dround pm3dGamma exp setgray} bind def\n");
      fprintf(fp, "} ifelse\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTB\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "210 4619 N\n");
      fprintf(fp, "210 168 L\n");
      fprintf(fp, "5802 0 V\n");
      fprintf(fp, "0 4451 V\n");
      fprintf(fp, "-5802 0 V\n");
      fprintf(fp, "Z stroke\n");
      fprintf(fp, "1.000 UP\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "3111 4829 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (Contact map distribution. accn:%s  (Res: %d", accn, nres);
      if(strcmp(chain, "-dummyval") != 0){
	    fprintf(fp, ", Chain: %s", chain);
      }
      if(nmr_mdl < 0){
	    fprintf(fp, ")");
      }else{
	    fprintf(fp, ", NMR MDL: %d)", nmr_mdl);
      }
      fprintf(fp, ")]\n");
      fprintf(fp, "] -46.7 MCshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "/vshift -46 def\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "%% Begin plot #1\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "%%%%%%%%BeginImage\n");
      fprintf(fp, "gsave 210 4619 N 210 168 L 6012 168 L 6012 4619 L Z clip\n");
      fprintf(fp, "InterpretLevel1 {\n");
      fprintf(fp, "  %%%% Construct a box instead of image\n");
      fprintf(fp, "  LTb\n");
      fprintf(fp, "  210 4619 M\n");
      fprintf(fp, "  5802 0 V\n");
      fprintf(fp, "  0 -4451 V\n");
      fprintf(fp, "  -5802 0 V\n");
      fprintf(fp, "  210 4619 L\n");
      fprintf(fp, "  40 -110 R\n");
      fprintf(fp, "  (PS level 2 image) Lshow\n");
      fprintf(fp, "  %% Read data but ignore it\n");
      int bytes = (nres/2)* (5*nres)/4;
      if (bytes > 65535) {
	    fprintf(fp, "  /imagebuf 65535 string def\n");
	    fprintf(fp, "  /imagebuf_rest %d string def\n", bytes % 65535);
	    fprintf(fp, "   1 1 %d { pop currentfile imagebuf readstring } for\n", bytes / 65535);
	    fprintf(fp, "  currentfile imagebuf_rest readstring\n");
      } else {
	    fprintf(fp, "  /imagebuf %d string def\n", bytes);
	    fprintf(fp, "  currentfile imagebuf readstring\n");
      }
      fprintf(fp, "} {\n");
      fprintf(fp, "gsave\n");
      fprintf(fp, "210 4619 translate\n");
      fprintf(fp, "5802 -4451 scale\n");
      fprintf(fp, "%%%%%%%%BeginPalette\n");
      fprintf(fp, "[ /Indexed\n");
      fprintf(fp, "  /DeviceRGB 8\n");
      fprintf(fp, "  <\n");
      fprintf(fp, "   252121 ff2bd5 a8a8a8 ffbc00 ff0000 0055ff 00ff00 ffffff\n");
      fprintf(fp, "   000000\n");
      fprintf(fp, "  >\n");
      fprintf(fp, "] setcolorspace\n");
      fprintf(fp, "%%%%%%%%EndPalette\n");
      fprintf(fp, "<<\n");
      fprintf(fp, "  /ImageType 1\n");
      fprintf(fp, "  /Width %d\n", nres);
      fprintf(fp, "  /Height %d\n", nres);
      fprintf(fp, "  /BitsPerComponent 4\n");
      fprintf(fp, "  /ImageMatrix [ %d 0 0 %d 0 0 ]\n", nres, nres);
      fprintf(fp, "  /Decode [ 0 15 ]\n");
      fprintf(fp, "  /DataSource currentfile /ASCII85Decode filter\n");
      fprintf(fp, "  /MultipleDataSources false\n");
      fprintf(fp, "  /Interpolate false\n");
      fprintf(fp, ">>\n");
      fprintf(fp, "image\n");
      fprintf(fp, "} ifelse\n");
}

void ps_heatmap_tail(FILE* fp)
{
      fprintf(fp, "InterpretLevel1 not {\n");
      fprintf(fp, "  grestore\n");
      fprintf(fp, "} if\n");
      fprintf(fp, "grestore\n");
      fprintf(fp, "%%%%%%%%EndImage\n");
      fprintf(fp, "%% End plot #1\n");
      fprintf(fp, "2.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTB\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "210 4619 N\n");
      fprintf(fp, "210 168 L\n");
      fprintf(fp, "5802 0 V\n");
      fprintf(fp, "0 4451 V\n");
      fprintf(fp, "-5802 0 V\n");
      fprintf(fp, "Z stroke\n");
      fprintf(fp, "stroke gsave	%%%% draw gray scale smooth box\n");
      fprintf(fp, "maxcolors 0 gt {/imax maxcolors def} {/imax 1024 def} ifelse\n");
      fprintf(fp, "6157 168 translate 290 4451 scale 0 setlinewidth\n");
      fprintf(fp, "/ystep 1 imax div def /y0 0 def /ii 0 def\n");
      fprintf(fp, "{ y0 g 0 y0 N 1 0 V 0 ystep V -1 0 f\n");
      fprintf(fp, "/y0 y0 ystep add def /ii ii 1 add def\n");
      fprintf(fp, "ii imax ge {exit} if } loop\n");
      fprintf(fp, "grestore 0 setgray\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6157 168 N\n");
      fprintf(fp, "290 0 V\n");
      fprintf(fp, "0 4451 V\n");
      fprintf(fp, "-290 0 V\n");
      fprintf(fp, "0 -4451 V\n");
      fprintf(fp, "Z stroke\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 415 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 415 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (BGCOL)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 909 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 909 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (ADJA)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 1404 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 1404 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (STACK)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 1898 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 1898 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (CROSS)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 2393 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 2393 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (NC-BP)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 2888 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 2888 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (WC-BP)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 3382 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 3382 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (CLOSE)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 3877 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 3877 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (PROX)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "6447 4371 M\n");
      fprintf(fp, "63 0 V\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "6594 4371 M\n");
      fprintf(fp, "[ [(Helvetica) 140.0 0.0 true true 0 (CHN)]\n");
      fprintf(fp, "] -46.7 MLshow\n");
      fprintf(fp, "/Helvetica findfont 140 scalefont setfont\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "1.000 UP\n");
      fprintf(fp, "1.000 UL\n");
      fprintf(fp, "LTb\n");
      fprintf(fp, "LCb setrgbcolor\n");
      fprintf(fp, "grestore %% colour palette end\n");
      fprintf(fp, "stroke\n");
      fprintf(fp, "grestore\n");
      fprintf(fp, "end\n");
      fprintf(fp, "showpage\n");
      fprintf(fp, "%%%%Trailer\n");
      fprintf(fp, "%%%%DocumentFonts: Helvetica\n");
      fprintf(fp, "%%%%Pages: 1\n");
}
void ps_heatmap_data(FILE* fp, int** mat, int nres)
{
	    int charcnt = 0;
	    const int maxchar=78;
	    unsigned char binary[4];
	    unsigned char data[8];
	    char asc85[6]; 
	    int cnt = 0;
	    int asclen = 0;

	    for(int i=nres-1; i>=0;  --i){
		  for(int j=0; j<nres; ++j){
			data[cnt] = (unsigned char) mat[i][j];
			cnt ++;
			if(cnt == 8){
			      cnt = 0;
			      tobinary(binary, data, 4);
			      toascii85(asc85, &asclen, binary, 4);
			      asc85[asclen] = '\0';
			      if(charcnt + asclen < maxchar){
				    fprintf(fp,"%s", asc85);
				    charcnt += asclen;

			      }else{
				    for(int k=0; k<asclen; ++k){
					  if(charcnt == maxchar){
						fprintf(fp, "\n");
						charcnt = 0;
					  }
					  fprintf(fp, "%c", asc85[k]);
					  charcnt ++;
				    }
			      }
			}
		  }
		  if(nres%2 == 1){
			data[cnt] = (unsigned char) 0;
			cnt ++;
			if(cnt == 8){
			      cnt = 0;
			      tobinary(binary, data, 4);
			      toascii85(asc85, &asclen, binary, 4);
			      asc85[asclen] = '\0';
			      if(charcnt + asclen < maxchar){
				    fprintf(fp,"%s", asc85);
				    charcnt += asclen;

			      }else{
				    for(int k=0; k<asclen; ++k){
					  if(charcnt == maxchar){
						fprintf(fp, "\n");
						charcnt = 0;
					  }
					  fprintf(fp, "%c", asc85[k]);
					  charcnt ++;
				    }
			      }
			}

		  }
	    }

	    if(cnt != 0){
		  for(int k=cnt; k<8; ++k){
			data[k] = 0;
		  }
		  tobinary(binary, data, 4);
		  toascii85(asc85, &asclen, binary, 4);
		  if(asclen ==1 && asc85[0] == 'z'){
			strcpy(asc85,"!!!!!");
		  }
		  asclen = (cnt+1)/2;
		  asc85[asclen+1] = '\0';

	    if(charcnt + asclen < maxchar){
		  fprintf(fp,"%s", asc85);

	    }else{
		  for(int k=0; k<asclen; ++k){
			if(charcnt == maxchar){
			      fprintf(fp, "\n");
			      charcnt = 0;
			}
			fprintf(fp, "%c", asc85[k]);
			charcnt ++;
		  }
	    }
	    }
	    fprintf(fp, "~>\n");
}
void ps_create_heatmap(const char* psfile, char* accn, int** mat, int nres, char* chain, char* nmrmodel)
{
      FILE* fp = fopen(psfile, "w");
      if(fp == NULL){    /* Exception Handling */ 
	    fprintf(stderr, "Error in function %s. (File: %s, Line %d)... Unable to create postsctipt file.\n", __func__, __FILE__, __LINE__);
	    exit(EXIT_FAILURE);
      }
      int nmrval;
      if(strcmp(nmrmodel, "-dummyval") ==0 ){
	    nmrval = -999;
      }else{
	    nmrval = atoi(nmrmodel);
      }

      ps_heatmap_head(fp, accn, nres, chain, nmrval);
      ps_heatmap_data(fp,mat,nres);
      ps_heatmap_tail(fp);
      fclose(fp);

}

#endif //__BIOLIB__POSTSCRIPT_H__
