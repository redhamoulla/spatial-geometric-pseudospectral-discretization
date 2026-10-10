#!/usr/bin/env python3
"""One command reproduces every committed numerical result (no network needed)."""
import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from gpps import MixedDiscretization,ProjectedGauss,integrate,transmission_line,shallow_water
from gpps.polynomials import gauss_rule
from gpps.spatial import tensor_diagnostics
from gpps.systems import PHSystem,interconnect
from gpps.time import exact_linear_flow


def table(out,name,rows):
    with (out/"tables"/f"{name}.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def physical_error(mesh,state,target):
    z,w=gauss_rule(60,mesh.interval)
    V=mesh.flow.values(z)
    exact=target(z)
    return float(np.sqrt(w@((V@state[:mesh.n]-exact[0])**2+(V@state[mesh.n:]-exact[1])**2)))


def rk4_reference(system,x0,T,steps,u=None):
    """Independent classical explicit RK4 reference, verified by doubling steps."""
    x=x0.copy();h=T/steps
    for n in range(steps):
        t=n*h
        k1=system.rhs(x,t,u);k2=system.rhs(x+h*k1/2,t+h/2,u)
        k3=system.rhs(x+h*k2/2,t+h/2,u);k4=system.rhs(x+h*k3,t+h,u)
        x+=h*(k1+2*k2+2*k3+k4)/6
    return x


def spatial_experiments(out,summary):
    rows=[]
    for family in ("gauss","lobatto","uniform","jittered"):
        for n in (2,4,8,12,16):
            m=MixedDiscretization(n,interval=(0,2),family=family)
            rows.append({"family":family,"n":n,**m.diagnostics()})
    table(out,"spatial_structure",rows)
    summary["spatial_dirac_max_n_le_12"]=max(r["dirac_relative"] for r in rows if r["n"]<=12)
    summary["spatial_dirac_max_n16"]=max(r["dirac_relative"] for r in rows if r["n"]==16)
    printed=np.array([.32063745754047,.96191238151097,1.60321563740732,2.24771134219675,
                      2.94830702957369,3.98374558826630,6.38606332778657,18.76179259703309])
    rows=[]
    for family in ("gauss","lobatto","uniform","jittered"):
        m=MixedDiscretization(8,interval=(0,2),family=family)
        system=transmission_line(m)
        ev=np.linalg.eigvals(system.J@system.hamiltonian.Q)
        omega=np.sort(ev.imag[ev.imag>0])
        for j,(om,ref) in enumerate(zip(omega,printed)):
            exact=(j+.5)*np.pi/(2*np.sqrt(6))
            rows.append({"family":family,"mode":j+1,"computed":om,"paper_table_4":ref,
                         "absolute_reproduction_error":abs(om-ref),"continuum":exact,
                         "relative_continuum_error":abs(om/exact-1),"max_real_part":np.max(np.abs(ev.real))})
    table(out,"paper_spectrum",rows)
    summary["paper_spectrum_max_abs_error"]=max(r["absolute_reproduction_error"] for r in rows)
    # Arbitrary coordinate changes must transport Hamiltonians, not just eigenvalues.
    a=MixedDiscretization(8,family="gauss");b=MixedDiscretization(8,family="jittered")
    T0=a.flow.values(b.z);T=np.block([[T0,np.zeros_like(T0)],[np.zeros_like(T0),T0]])
    sa,sb=transmission_line(a),transmission_line(b)
    invT=np.linalg.solve(T,np.eye(16))
    summary["coordinate_J_relative"]=float(np.linalg.norm(sb.J-T@sa.J@T.T)/np.linalg.norm(sb.J))
    summary["coordinate_Q_relative"]=float(np.linalg.norm(sb.hamiltonian.Q-invT.T@sa.hamiltonian.Q@invT)/np.linalg.norm(sb.hamiltonian.Q))
    rows=[]
    for nx in (2,4,8):
        for nt in (1,2,3):
            rows.append({"space_degree":nx,"time_degree":nt,**tensor_diagnostics(MixedDiscretization(nx),MixedDiscretization(nt))})
    table(out,"space_time_geometry",rows)
    summary["tensor_stokes_max"]=max(r["stokes_relative"] for r in rows)
    summary["tensor_complex_max"]=max(r["complex_relative"] for r in rows)


def convergence_experiments(out,summary):
    m=MixedDiscretization(8,interval=(0,2));sys=transmission_line(m)
    k=3*np.pi/4;x=np.r_[np.sin(k*m.z),np.zeros(8)]
    T=2.0;target=exact_linear_flow(sys,x,T);rows=[]
    for degree in (1,2):
        previous=None
        for steps in (10,20,40,80,160):
            path,ledger=integrate(sys,x,T/steps,steps,ProjectedGauss(degree,degree))
            error=float(np.sqrt((path[-1]-target)@sys.hamiltonian.Q@(path[-1]-target)))
            order=None if previous is None else float(np.log2(previous/error))
            rows.append({"degree":degree,"dt":T/steps,"steps":steps,"error":error,
                         "observed_order":order,"energy_defect":max(abs(r["balance"]) for r in ledger)})
            previous=error
        summary[f"linear_order_s{degree}"]=order
    table(out,"time_convergence",rows)
    rows=[]
    for n in (2,3,4,6,8,10,12):
        mesh=MixedDiscretization(n,interval=(0,2));system=transmission_line(mesh)
        xx=np.r_[np.sin(k*mesh.z),np.zeros(n)]
        exact=exact_linear_flow(system,xx,T)
        omega=k/np.sqrt(6)
        target=lambda z:(np.sin(k*z)*np.cos(omega*T),-np.sqrt(2/3)*np.cos(k*z)*np.sin(omega*T))
        path,_=integrate(system,xx,.025,80,ProjectedGauss(2,2))
        rows.append({"n":n,"spatial_error_exact_time":physical_error(mesh,exact,target),
                     "space_time_error_dt_0025":physical_error(mesh,path[-1],target)})
    table(out,"space_convergence",rows)
    summary["space_error_n12"]=rows[-1]["spatial_error_exact_time"]
    summary["space_time_error_n12"]=rows[-1]["space_time_error_dt_0025"]


def nonlinear_experiments(out,summary):
    m=MixedDiscretization(6);sys=transmission_line(m,1,1,quartic=1)
    x=np.r_[.65*np.sin(np.pi*m.z/2),.2*np.cos(np.pi*m.z/2)]
    all_rows=[]
    for name,s,k in (("midpoint",1,1),("AVF",1,2),("Gauss2",2,2),("projected2",2,4)):
        _,rows=integrate(sys,x,.03,1000,ProjectedGauss(s,k))
        h0=rows[0]["energy"]
        for r in rows:
            r["method"]=name;r["relative_energy_error"]=(r["energy"]-h0)/h0
        table(out,f"nonlinear_{name}",rows)
        summary[f"nonlinear_energy_{name}"]=max(abs(r["relative_energy_error"]) for r in rows)
        all_rows.extend(rows)
    # An independent trajectory reference distinguishes accuracy from conservation.
    T=.6;ref1=rk4_reference(sys,x,T,6000);ref2=rk4_reference(sys,x,T,12000)
    summary["nonlinear_reference_refinement"]=float(np.linalg.norm(ref1-ref2))
    rows=[]
    for s,k in ((1,2),(2,4)):
        previous=None
        for steps in (20,40,80,160):
            path,_=integrate(sys,x,T/steps,steps,ProjectedGauss(s,k))
            error=float(np.linalg.norm(path[-1]-ref2))
            order=None if previous is None else float(np.log2(previous/error))
            rows.append({"degree":s,"dt":T/steps,"error":error,"observed_order":order})
            previous=error
        summary[f"nonlinear_order_s{s}"]=order
    table(out,"nonlinear_convergence",rows)
    damped=transmission_line(m,1,1,resistance=.15,quartic=1)
    u=lambda t:np.array([.2*np.sin(2*t)*np.exp(-.15*t),.12*np.sin(3*t)])
    _,rows=integrate(damped,x,.02,500,ProjectedGauss(2,4),u)
    table(out,"forced_damped",rows)
    summary["forced_damped_balance"]=max(abs(r["balance"]) for r in rows)
    summary["forced_damped_step_balance"]=max(abs(r["step_balance"]) for r in rows)
    summary["forced_damped_total_dissipation"]=rows[-1]["dissipation"]
    summary["forced_damped_min_dissipation_increment"]=min(rows[i]["dissipation"]-rows[i-1]["dissipation"] for i in range(1,len(rows)))


def shallow_water_experiment(out,summary):
    m=MixedDiscretization(12);sys=shallow_water(m)
    x=np.r_[1+.02*np.sin(np.pi*m.z/2),np.zeros(12)]
    u=lambda t:np.array([1.,0.])
    path,rows=integrate(sys,x,.002,300,ProjectedGauss(2,3),u)
    z=np.linspace(0,1,401);V=m.flow.values(z)
    depths=path[:,:12]@V.T
    table(out,"shallow_water_ledger",rows)
    profiles=[]
    for i in (0,100,200,300):
        for zz,hh,vv in zip(z,depths[i],V@path[i,12:]):
            profiles.append({"time":i*.002,"z":zz,"depth":hh,"velocity":vv})
    table(out,"shallow_water_profiles",profiles)
    summary["shallow_water_balance"]=max(abs(r["balance"]) for r in rows)
    summary["shallow_water_minimum_sampled_depth"]=float(depths.min())
    rest=np.r_[np.ones(12),np.zeros(12)]
    summary["shallow_water_rest_residual"]=float(np.linalg.norm(sys.rhs(rest,0,u)))
    # Short smooth interval: compare two space resolutions, not a shock benchmark.
    fine=MixedDiscretization(18);fsys=shallow_water(fine)
    fx=np.r_[1+.02*np.sin(np.pi*fine.z/2),np.zeros(18)]
    fp,_=integrate(fsys,fx,.001,600,ProjectedGauss(2,3),u)
    summary["shallow_water_final_refinement_difference"]=float(np.sqrt(np.mean((V@path[-1,:12]-fine.flow.values(z)@fp[-1,:18])**2)))


def network_experiment(out,summary):
    ma,mb=MixedDiscretization(6),MixedDiscretization(5)
    a,b=transmission_line(ma,1,1),transmission_line(mb,1.4,.8)
    K=np.zeros((4,4));K[1,2]=1;K[2,1]=-1
    sys,ports=interconnect([a,b],K)
    x=np.r_[.4*np.sin(np.pi*ma.z/2),np.zeros(6),np.zeros(10)]
    power=[]
    def observe(i,step):
        for e in step.stage_efforts:
            u,y=ports(e,np.zeros(4))
            power.append(abs(u[1]*y[1]+u[2]*y[2]))
    path,rows=integrate(sys,x,.01,600,ProjectedGauss(2,2),observer=observe)
    for x,r in zip(path,rows):
        r["energy_line_a"]=a.hamiltonian.value(x[:12])
        r["energy_line_b"]=b.hamiltonian.value(x[12:])
    table(out,"interconnection",rows)
    summary["network_balance"]=max(abs(r["balance"]) for r in rows)
    summary["network_internal_power_cancellation"]=max(power)
    summary["network_max_energy_line_b"]=max(r["energy_line_b"] for r in rows)


def oscillator_experiment(out,summary):
    class Energy:
        degree=4
        def value(self,x):return .5*np.sum(x*x)+.25*x[0]**4
        def gradient(self,x):return np.array([x[0]+x[0]**3,x[1]])
        def hessian(self,x):return np.diag([1+3*x[0]**2,1])
    sys=PHSystem(np.array([[0.,1.],[-1.,0.]]),np.zeros((2,2)),np.zeros((2,0)),np.zeros((0,0)),Energy())
    rows=[]
    for name,s,k in (("midpoint",1,1),("AVF",1,2),("Gauss2",2,2),("projected2",2,4)):
        result=ProjectedGauss(s,k).step(sys,np.array([.8,.7]),0,.2,tangent=True)
        rows.append({"method":name,"energy_defect":abs(result.balance_residual),
                     "symplectic_defect":float(np.linalg.norm(result.tangent.T@sys.J@result.tangent-sys.J)),
                     "chain_defect":abs(result.chain_defect)})
    table(out,"symplectic_vs_energy",rows)
    summary["midpoint_symplectic_defect"]=rows[0]["symplectic_defect"]
    summary["avf_symplectic_defect"]=rows[1]["symplectic_defect"]


def make_figures(out):
    from figures import plot
    def read(name):
        with (out/"tables"/f"{name}.csv").open() as f:return list(csv.DictReader(f))
    def points(rows,x,y):return [(float(r[x]),float(r[y])) for r in rows]
    rows=read("time_convergence")
    series=[]
    for s in (1,2):
        series.append((f"degree {s}",[(np.log10(float(r["dt"])),np.log10(float(r["error"]))) for r in rows if int(r["degree"])==s]))
    plot(out/"figures"/"time_convergence.svg","Temporal convergence","Closed transmission line; independent matrix-exponential reference",series,"log10(time step)","log10(energy-norm error)")
    series=[]
    for name in ("midpoint","AVF","Gauss2","projected2"):
        rows=read(f"nonlinear_{name}")
        series.append((name,[(float(r["time"]),np.log10(max(abs(float(r["relative_energy_error"])),1e-16))) for r in rows[1:]]))
    plot(out/"figures"/"nonlinear_energy.svg","Nonlinear energy preservation","Quartic distributed Hamiltonian; dt = 0.03, T = 30",series,"time","log10(relative energy error); display floor 1e-16")
    rows=read("forced_damped");h0=float(rows[0]["energy"])
    series=[("H(t) - H(0)",[(float(r["time"]),float(r["energy"])-h0) for r in rows]),
            ("boundary work",points(rows,"time","work")),("minus dissipation",[(float(r["time"]),-float(r["dissipation"])) for r in rows]),
            ("work minus dissipation",[(float(r["time"]),float(r["work"])-float(r["dissipation"])) for r in rows])]
    plot(out/"figures"/"forced_balance.svg","Open-system energy ledger","Nonlinear line, two boundary inputs and physical resistance",series,"time","energy / work")
    rows=read("shallow_water_profiles")
    series=[(f"t = {t:.1f}",points([r for r in rows if abs(float(r["time"])-t)<1e-10],"z","depth")) for t in (0.,.2,.4,.6)]
    plot(out/"figures"/"shallow_water.svg","Nonlinear shallow-water dynamics","Smooth small-amplitude flow; 12 spatial states per field",series,"position","water depth")
    rows=read("interconnection")
    plot(out/"figures"/"interconnection.svg","Energy exchange between two lines","Power-conserving interface, including boundary feedthrough",[("line A",points(rows,"time","energy_line_a")),("line B",points(rows,"time","energy_line_b")),("total",points(rows,"time","energy"))],"time","energy")
    rows=read("space_convergence")
    plot(out/"figures"/"space_convergence.svg","Spatial and space-time accuracy","Exact standing-wave reference; fixed dt = 0.025 for the full simulation",[("exact time evolution",[(float(r["n"]),np.log10(float(r["spatial_error_exact_time"]))) for r in rows]),("degree 2 time scheme",[(float(r["n"]),np.log10(float(r["space_time_error_dt_0025"]))) for r in rows])],"number of spatial states per field","log10(L2 field error)")


def validate(s):
    checks={
        "spectrum reproduced":s["paper_spectrum_max_abs_error"]<1e-8,
        "Dirac structure":s["spatial_dirac_max_n_le_12"]<1e-8,
        "space-time Stokes":s["tensor_stokes_max"]<1e-11,
        "exterior complex":s["tensor_complex_max"]<1e-13,
        "linear order 2":s["linear_order_s1"]>1.9,
        "linear order 4":s["linear_order_s2"]>3.8,
        "nonlinear order 2":s["nonlinear_order_s1"]>1.8,
        "nonlinear order 4":s["nonlinear_order_s2"]>3.5,
        "independent nonlinear reference":s["nonlinear_reference_refinement"]<1e-10,
        "nonlinear energy":s["nonlinear_energy_projected2"]<1e-9,
        "AVF energy":s["nonlinear_energy_AVF"]<1e-9,
        "negative energy control":s["nonlinear_energy_midpoint"]>1e-5,
        "forced balance":s["forced_damped_balance"]<1e-9,
        "nonnegative dissipation":s["forced_damped_min_dissipation_increment"]>=0,
        "shallow water balance":s["shallow_water_balance"]<1e-9,
        "positive sampled depth":s["shallow_water_minimum_sampled_depth"]>.9,
        "rest state":s["shallow_water_rest_residual"]<1e-9,
        "network balance":s["network_balance"]<1e-9,
        "internal power cancellation":s["network_internal_power_cancellation"]<1e-10,
        "nontrivial energy transfer":s["network_max_energy_line_b"]>.01,
        "symplectic control":s["midpoint_symplectic_defect"]<1e-12,
        "AVF not symplectic":s["avf_symplectic_defect"]>1e-6,
    }
    return {key: bool(value) for key,value in checks.items()}


def main():
    p=argparse.ArgumentParser();p.add_argument("--output",type=Path,default=ROOT/"results")
    p.add_argument("--no-figures",action="store_true");p.add_argument("--check",action="store_true")
    args=p.parse_args();out=args.output
    (out/"tables").mkdir(parents=True,exist_ok=True);(out/"figures").mkdir(exist_ok=True)
    summary={};t0=time.perf_counter()
    for task in (spatial_experiments,convergence_experiments,nonlinear_experiments,shallow_water_experiment,network_experiment,oscillator_experiment):
        print(task.__name__,flush=True);task(out,summary)
    checks=validate(summary)
    summary["checks"]=checks
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    manifest={"python":platform.python_version(),"numpy":np.__version__,"platform":platform.system()+" "+platform.machine(),
              "elapsed_seconds":time.perf_counter()-t0,"seed":20261010,"sources":{}}
    for file in sorted(list((ROOT/"src").rglob("*.py"))+list((ROOT/"experiments").glob("*.py"))):
        manifest["sources"][str(file.relative_to(ROOT))]=hashlib.sha256(file.read_bytes()).hexdigest()
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    if not args.no_figures:make_figures(out)
    print(json.dumps(summary,indent=2),flush=True)
    if args.check and not all(checks.values()):
        raise SystemExit("Failed validation gates: "+", ".join(k for k,v in checks.items() if not v))


if __name__=="__main__":main()
