F32 eps2

rij = EPJ.pos - EPI.pos
r2  = rij * rij + eps2
r_inv  = rsqrt(r2)
tmp    = 3.0f - r2*(r_inv*r_inv)
r_inv *= (tmp * 0.5f)

r2_inv  = r_inv  * r_inv
r3_inv  = r2_inv * r_inv
r4_inv  = r2_inv * r2_inv
r5_inv  = r2_inv * r3_inv

//tr  = qj_xxloc + qj_yyloc + qj_zzloc
//qxx = 3.0f * qj_xxloc - tr
//qyy = 3.0f * qj_yyloc - tr
//qzz = 3.0f * qj_zzloc - tr
//qxy = 3.0f * qj_xyloc
//qyz = 3.0f * qj_yzloc
//qzx = 3.0f * qj_zxloc
//mtr = -(eps2 * tr)

tr  = EPJ.quad_xx + EPJ.quad_yy + EPJ.quad_zz
qxx = 3.0f * EPJ.quad_xx - tr
qyy = 3.0f * EPJ.quad_yy - tr
qzz = 3.0f * EPJ.quad_zz - tr
qxy = 3.0f * EPJ.quad_xy
qyz = 3.0f * EPJ.quad_yz
qzx = 3.0f * EPJ.quad_xz
mtr = -(eps2 * tr)

qr.x = qxx*rij.x + qxy*rij.y + qzx*rij.z
qr.y = qyy*rij.y + qyz*rij.z + qxy*rij.x
qr.z = qzz*rij.z + qzx*rij.x + qyz*rij.y
rqr  = mtr + qr * rij
rqr_r4_inv = rqr * r4_inv

meff  =  EPJ.mass + 0.5f * rqr_r4_inv
meff3 = (EPJ.mass + 2.5f * rqr_r4_inv) * r3_inv

FORCE.pot -= meff * r_inv
FORCE.acc -= r5_inv * qr + meff3 * rij


