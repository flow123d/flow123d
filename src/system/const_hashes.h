/*
 * const_hashes.h
 *
 *  Created on: Nov 5, 2012
 *      Author: jb
 */

#ifndef CONST_HASHES_H_
#define CONST_HASHES_H_


// maximal string length that influence hash - the depth of the hashing
#define CONSTHASH_DEPTH 64

// randomly generated constants.  The bottom half has to be FFFF or
// else the entire hash loses some strength
static const size_t CONSTHASH_CONSTANTS[CONSTHASH_DEPTH+1] =
{
    0x6157FFFF, 0x5387FFFF, 0x8ECBFFFF, 0xB3DBFFFF, 0x1AFDFFFF, 0xD1FDFFFF, 0x19B3FFFF, 0xE6C7FFFF,
    0x53BDFFFF, 0xCDAFFFFF, 0xE543FFFF, 0x369DFFFF, 0x8135FFFF, 0x50E1FFFF, 0x115BFFFF, 0x07D1FFFF,
    0x9AA1FFFF, 0x4D87FFFF, 0x442BFFFF, 0xEAA5FFFF, 0xAEDBFFFF, 0xB6A5FFFF, 0xAFE9FFFF, 0xE895FFFF,
    0x4E05FFFF, 0xF8BFFFFF, 0xCB5DFFFF, 0x2F85FFFF, 0xF1DFFFFF, 0xD5ADFFFF, 0x438DFFFF, 0x6073FFFF,
    0xA99FFFFF, 0x2E0BFFFF, 0xF729FFFF, 0x5D01FFFF, 0x1ACDFFFF, 0xFAD1FFFF, 0xD86BFFFF, 0xE909FFFF,
    0xD3BDFFFF, 0xF35BFFFF, 0xD53DFFFF, 0x4DC1FFFF, 0x6897FFFF, 0x6E4DFFFF, 0x305BFFFF, 0x6F0DFFFF,
    0x33C9FFFF, 0xC955FFFF, 0xC1EDFFFF, 0x48D5FFFF, 0x0CF5FFFF, 0x356BFFFF, 0x5F65FFFF, 0x71C1FFFF,
    0x3F13FFFF, 0x489DFFFF, 0xEBA3FFFF, 0x340DFFFF, 0xF537FFFF, 0xD5E7FFFF, 0x6D27FFFF, 0x89D7FFFF,
    0xA93FFFFF,
};

// multiplication constants, this allows an abstract use
// of the string length
static const size_t CONSTHASH_MULTS[CONSTHASH_DEPTH+1] =
{
    33,  34,  35,  36,  37,  38,  39,  40,
    41,  42,  43,  44,  45,  46,  47,  48,
    49,  50,  51,  52,  53,  54,  55,  56,
    57,  58,  59,  60,  61,  62,  63,  64,
    65,  66,  67,  68,  69,  70,  71,  72,
    73,  74,  75,  76,  77,  78,  79,  80,
    81,  82,  83,  84,  85,  86,  87,  88,
    89,  90,  91,  92,  93,  94,  95,  96,
    97,
};

#define CONSTHASH_RECURSE_00(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[0] :  CONSTHASH_RECURSE_01(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_01(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[1] :  CONSTHASH_RECURSE_02(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_02(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[2] :  CONSTHASH_RECURSE_03(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_03(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[3] :  CONSTHASH_RECURSE_04(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_04(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[4] :  CONSTHASH_RECURSE_05(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_05(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[5] :  CONSTHASH_RECURSE_06(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_06(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[6] :  CONSTHASH_RECURSE_07(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_07(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[7] :  CONSTHASH_RECURSE_08(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_08(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[8] :  CONSTHASH_RECURSE_09(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_09(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[9] :  CONSTHASH_RECURSE_10(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_10(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[10] :  CONSTHASH_RECURSE_11(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_11(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[11] :  CONSTHASH_RECURSE_12(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_12(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[12] :  CONSTHASH_RECURSE_13(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_13(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[13] :  CONSTHASH_RECURSE_14(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_14(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[14] :  CONSTHASH_RECURSE_15(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_15(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[15] :  CONSTHASH_RECURSE_16(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_16(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[16] :  CONSTHASH_RECURSE_17(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_17(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[17] :  CONSTHASH_RECURSE_18(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_18(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[18] :  CONSTHASH_RECURSE_19(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_19(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[19] :  CONSTHASH_RECURSE_20(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_20(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[20] :  CONSTHASH_RECURSE_21(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_21(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[21] :  CONSTHASH_RECURSE_22(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_22(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[22] :  CONSTHASH_RECURSE_23(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_23(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[23] :  CONSTHASH_RECURSE_24(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_24(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[24] :  CONSTHASH_RECURSE_25(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_25(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[25] :  CONSTHASH_RECURSE_26(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_26(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[26] :  CONSTHASH_RECURSE_27(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_27(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[27] :  CONSTHASH_RECURSE_28(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_28(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[28] :  CONSTHASH_RECURSE_29(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_29(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[29] :  CONSTHASH_RECURSE_30(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_30(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[30] :  CONSTHASH_RECURSE_31(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_31(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[31] :  CONSTHASH_RECURSE_32(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_32(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[32] :  CONSTHASH_RECURSE_33(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_33(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[33] :  CONSTHASH_RECURSE_34(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_34(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[34] :  CONSTHASH_RECURSE_35(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_35(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[35] :  CONSTHASH_RECURSE_36(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_36(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[36] :  CONSTHASH_RECURSE_37(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_37(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[37] :  CONSTHASH_RECURSE_38(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_38(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[38] :  CONSTHASH_RECURSE_39(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_39(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[39] :  CONSTHASH_RECURSE_40(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_40(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[40] :  CONSTHASH_RECURSE_41(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_41(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[41] :  CONSTHASH_RECURSE_42(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_42(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[42] :  CONSTHASH_RECURSE_43(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_43(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[43] :  CONSTHASH_RECURSE_44(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_44(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[44] :  CONSTHASH_RECURSE_45(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_45(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[45] :  CONSTHASH_RECURSE_46(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_46(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[46] :  CONSTHASH_RECURSE_47(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_47(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[47] :  CONSTHASH_RECURSE_48(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_48(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[48] :  CONSTHASH_RECURSE_49(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_49(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[49] :  CONSTHASH_RECURSE_50(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_50(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[50] :  CONSTHASH_RECURSE_51(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_51(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[51] :  CONSTHASH_RECURSE_52(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_52(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[52] :  CONSTHASH_RECURSE_53(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_53(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[53] :  CONSTHASH_RECURSE_54(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_54(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[54] :  CONSTHASH_RECURSE_55(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_55(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[55] :  CONSTHASH_RECURSE_56(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_56(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[56] :  CONSTHASH_RECURSE_57(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_57(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[57] :  CONSTHASH_RECURSE_58(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_58(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[58] :  CONSTHASH_RECURSE_59(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_59(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[59] :  CONSTHASH_RECURSE_60(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_60(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[60] :  CONSTHASH_RECURSE_61(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_61(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[61] :  CONSTHASH_RECURSE_62(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_62(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[62] :  CONSTHASH_RECURSE_63(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_63(string, value) CONSTHASH_FUNCTION((*(string+1) == 0 ? CONSTHASH_CONSTANTS[63] :  CONSTHASH_RECURSE_64(string+1, *(string+1))), value)
#define CONSTHASH_RECURSE_64(string, value) CONSTHASH_CONSTANTS[64]

// The following is the function used for hashing
// Do NOT use NEXTHASH more than once, it will cause
// N-Squared expansion and make compilation very slow
// If not impossible

#define CONSTHASH_FUNCTION(next, value) ( (next) * 101 + (unsigned int) (value) );

// finally the macro used to generate the hash
#define CONSTHASH(string) CONSTHASH_RECURSE_00(string, *string)



#endif /* CONST_HASHES_H_ */
