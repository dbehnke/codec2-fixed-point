#ifndef CODEC2_MOD_QUANTISE_H
#define CODEC2_MOD_QUANTISE_H

#include <stdint.h>

void pack(uint8_t *bitArray, unsigned int *bitIndex,
          unsigned int field, unsigned int fieldWidth);

int unpack(const uint8_t *bitArray, unsigned int *bitIndex,
           unsigned int fieldWidth);

int decode_Wo_q(int index);
uint32_t wo_q_to_angle(int wq);
int encode_energy_log2(int32_t log2E_q16);
int32_t decode_energy_log2_q16(int index);
void encode_lspds_q16(int *indexes, const int32_t *lsp_hz_q16);
void decode_lspds_hz2(int16_t *hz2, const int *indexes);

#endif /* CODEC2_MOD_QUANTISE_H */
