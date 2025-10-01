#include "continuityErrs.H"
#include "adjustPhi.H"

// Adjust face flux to be consistent with p & U BCs
adjustPhi(phi, U, p);
// If your solver uses p_rgh instead of p, use:
// adjustPhi(phi, U, p_rgh);
