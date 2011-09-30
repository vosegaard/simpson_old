#ifndef __CMBLOCK_H
#define __CMBLOCK_H

int m_getnqlist();
void m_destroyQlist();
void m_makeQlist(mv_complx* Q);
void m_simtransadjQ(mv_complx * L,mv_complx* U);


#endif
