"""Local exact outcome-controlled CCZ through an ordinary measured AND.

Upstream's authoring API only permits outcome conditions on Clifford phases.
This emits one normal AND, a conditioned CZ and its measurement eraser, all
of which the original builder and native host already support and charge.
No upstream file, operation accounting, translator or verifier is changed.
"""
def install(m):
 Base=m.Builder
 class Builder(Base):
  def ccz(self,a,b,c,cond=None):
   if cond is None:return super().ccz(a,b,c)
   unique=list(dict.fromkeys((a,b,c)))
   if len(unique)==1:return self.z(unique[0],cond=cond)
   if len(unique)==2:return self.cz(*unique,cond=cond)
   work=self.qubit();self.ccx(a,b,work);self.cz(work,c,cond=cond)
   measurement=self.hmr(work);self.cz(a,b,cond=measurement);self.recycle(work)
 m.Builder=Builder

