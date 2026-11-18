#include "codec2_internal.h"
#include "interp.h"

/* Decoder interpolation of the 10 ms subframe between two 20 ms parameter sets, in integers.
 * Same decisions and values as the float interp_Wo / interp_energy / interpolate_lsp. */

int interp_Wo_q(int *voiced, int prev_voiced, int prev_wq, int next_voiced, int next_wq)
{
    /* trap corner case where voicing est is probably wrong */
    if (*voiced && !prev_voiced && !next_voiced)
        *voiced = 0;

    /* Wo depends on voicing of this and adjacent frames */
    if (!*voiced)
        return WO_Q_MIN;
    if (prev_voiced && next_voiced)
        return (prev_wq + next_wq) >> 1; /* both even (WO_Q_MIN, WO_Q_STEP even): exact */
    return prev_voiced ? prev_wq : next_wq;
}

int32_t interp_energy_log2_q16(int32_t prev_log2e, int32_t next_log2e)
{
    /* sqrt(prev_e next_e) in the log domain; |log2 E| < 2^20 */
    return (prev_log2e + next_log2e) >> 1;
}

void interpolate_lsp_hz2(int16_t *interp, const int16_t *restrict prev, const int16_t *restrict next)
{
    /* decoded LSPs and the initial state are even in half-Hz, so the mean is exact */
    for (int i = 0; i < LPC_ORD; i++)
        interp[i] = (int16_t)((prev[i] + next[i]) >> 1);
}
