#****************************************************************************                               
#                     Molecular Dynamics Potentials (MDP)
#                            CESMIX-MIT Project  
# 
#  Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
#****************************************************************************

using Libdl

global snaplib

function loadsnap(mdppath::String)    
    libexist = 0
    if Sys.isapple()
        libpath = mdppath * "src/Potential/cpuSnap.dylib";
    elseif Sys.isunix()
        libpath = mdppath * "src/Potential/cpuSnap.so";
    elseif Sys.iswindows()
        libpath = mdppath * "src/Potential/cpuSnap.dll";
    end
    if isfile(libpath)
        libexist = 1;
    end        
    if libexist==0            
        cdir = pwd();
        cd(mdppath * "src/Potential")        
        cstr = `g++ -std=c++11 -O3 -Wall -Wextra -pedantic -c -fPIC cpuSnap.cpp -o cpusnap.o`                
        run(cstr)
        if Sys.isapple()            
            cstr = `g++ -shared cpusnap.o -o cpuSnap.dylib`            
        elseif Sys.isunix()
            cstr = `g++ -shared cpusnap.o -o cpuSnap.so`     
        elseif Sys.iswindows()
            cstr = `g++ -shared cpusnap.o -o cpuSnap.dll`     
        end 
        run(cstr)
        cd(cdir)                   
    end
    global snaplib = Libdl.dlopen(libpath)
end

function closesnap()
    Libdl.dlclose(snaplib)
end

function bnum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxbcount), Cint, (Cint,), twojmax)
    return count
end

function unum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxucount), Cint, (Cint,), twojmax)
    return count
end

function znum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxzcount), Cint, (Cint,), twojmax)
    return count
end

function cgnum(twojmax::Int32)
    count = ccall(Libdl.dlsym(snaplib, :idxcgcount), Cint, (Cint,), twojmax)
    return count
end

function factable()

    facto = [
        1,
        1,
        2,
        6,
        24,
        120,
        720,
        5040,
        40320,
        362880,
        3628800,
        39916800,
        479001600,
        6227020800,
        87178291200,
        1307674368000,
        20922789888000,
        355687428096000,
        6.402373705728e+15,
        1.21645100408832e+17,
        2.43290200817664e+18,
        5.10909421717094e+19,
        1.12400072777761e+21,
        2.5852016738885e+22,
        6.20448401733239e+23,
        1.5511210043331e+25,
        4.03291461126606e+26,
        1.08888694504184e+28,
        3.04888344611714e+29,
        8.8417619937397e+30,
        2.65252859812191e+32,
        8.22283865417792e+33,
        2.63130836933694e+35,
        8.68331761881189e+36,
        2.95232799039604e+38,
        1.03331479663861e+40,
        3.71993326789901e+41,
        1.37637530912263e+43,
        5.23022617466601e+44,
        2.03978820811974e+46,
        8.15915283247898e+47,
        3.34525266131638e+49,
        1.40500611775288e+51,
        6.04152630633738e+52,
        2.65827157478845e+54,
        1.1962222086548e+56,
        5.50262215981209e+57,
        2.58623241511168e+59,
        1.24139155925361e+61,
        6.08281864034268e+62,
        3.04140932017134e+64,
        1.55111875328738e+66,
        8.06581751709439e+67,
        4.27488328406003e+69,
        2.30843697339241e+71,
        1.26964033536583e+73,
        7.10998587804863e+74,
        4.05269195048772e+76,
        2.35056133128288e+78,
        1.3868311854569e+80,
        8.32098711274139e+81,
        5.07580213877225e+83,
        3.14699732603879e+85,
        1.98260831540444e+87,
        1.26886932185884e+89,
        8.24765059208247e+90,
        5.44344939077443e+92,
        3.64711109181887e+94,
        2.48003554243683e+96,
        1.71122452428141e+98,
        1.19785716699699e+100,
        8.50478588567862e+101,
        6.12344583768861e+103,
        4.47011546151268e+105,
        3.30788544151939e+107,
        2.48091408113954e+109,
        1.88549470166605e+111,
        1.45183092028286e+113,
        1.13242811782063e+115,
        8.94618213078297e+116,
        7.15694570462638e+118,
        5.79712602074737e+120,
        4.75364333701284e+122,
        3.94552396972066e+124,
        3.31424013456535e+126,
        2.81710411438055e+128,
        2.42270953836727e+130,
        2.10775729837953e+132,
        1.85482642257398e+134,
        1.65079551609085e+136,
        1.48571596448176e+138,
        1.3520015276784e+140,
        1.24384140546413e+142,
        1.15677250708164e+144,
        1.08736615665674e+146,
        1.03299784882391e+148,
        9.91677934870949e+149,
        9.61927596824821e+151,
        9.42689044888324e+153,
        9.33262154439441e+155,
        9.33262154439441e+157,
        9.42594775983835e+159,
        9.61446671503512e+161,
        9.90290071648618e+163,
        1.02990167451456e+166,
        1.08139675824029e+168,
        1.14628056373471e+170,
        1.22652020319614e+172,
        1.32464181945183e+174,
        1.44385958320249e+176,
        1.58824554152274e+178,
        1.76295255109024e+180,
        1.97450685722107e+182,
        2.23119274865981e+184,
        2.54355973347219e+186,
        2.92509369349301e+188,
        3.3931086844519e+190,
        3.96993716080872e+192,
        4.68452584975429e+194,
        5.5745857612076e+196,
        6.68950291344912e+198,
        8.09429852527344e+200,
        9.8750442008336e+202,
        1.21463043670253e+205,
        1.50614174151114e+207,
        1.88267717688893e+209,
        2.37217324288005e+211,
        3.01266001845766e+213,
        3.8562048236258e+215,
        4.97450422247729e+217,
        6.46685548922047e+219,
        8.47158069087882e+221,
        1.118248651196e+224,
        1.48727070609069e+226,
        1.99294274616152e+228,
        2.69047270731805e+230,
        3.65904288195255e+232,
        5.01288874827499e+234,
        6.91778647261949e+236,
        9.61572319694109e+238,
        1.34620124757175e+241,
        1.89814375907617e+243,
        2.69536413788816e+245,
        3.85437071718007e+247,
        5.5502938327393e+249,
        8.04792605747199e+251,
        1.17499720439091e+254,
        1.72724589045464e+256,
        2.55632391787286e+258,
        3.80892263763057e+260,
        5.71338395644585e+262,
        8.62720977423323e+264,
        1.31133588568345e+267,
        2.00634390509568e+269,
        3.08976961384735e+271,
        4.78914290146339e+273,
        7.47106292628289e+275,
        1.17295687942641e+278,
        1.85327186949373e+280,
        2.94670227249504e+282,
        4.71472363599206e+284,
        7.59070505394721e+286,
        1.22969421873945e+289,
        2.0044015765453e+291,
        3.28721858553429e+293,
        5.42391066613159e+295,
        9.00369170577843e+297,
        1.503616514865e+300,
    ]
    
    facto = Float64.(Array(facto))
    
    return facto
    
end

function initsnap(snap::SNAPparams, snapcoeff=nothing)    
    param = zeros(11)
    param[1] = length(snap.species)
    param[2] = snap.twojmax  
    param[3] = snap.rcutfac 
    param[4] = snap.rfac0
    param[5] = snap.rmin0 
    param[6] = snap.bzeroflag
    param[7] = snap.switchflag
    param[8] = snap.quadraticflag
    param[9] = snap.chemflag
    param[10] = snap.bnormflag
    param[11] = snap.wselfallflag

    return initsna(param, snap.elemradius, snap.elemweight, snapcoeff)
end

function initsnap(snap::SNAPparams, elemradius, elemweight, snapcoeff=nothing)

    param = zeros(11)
    param[1] = length(snap.species)
    param[2] = snap.twojmax  
    param[3] = snap.rcutfac 
    param[4] = snap.rfac0
    param[5] = snap.rmin0 
    param[6] = snap.bzeroflag
    param[7] = snap.switchflag
    param[8] = snap.quadraticflag
    param[9] = snap.chemflag
    param[10] = snap.bnormflag
    param[11] = snap.wselfallflag

    return initsna(param, elemradius, elemweight, snapcoeff)
end

function initsna(param, elemradius, elemweight, snapcoeff=nothing)

    ntypes = Int32(param[1]);  
    twojmax = Int32(param[2]);  
    rcutfac = param[3];
    rfac0 = param[4];
    rmin0 = param[5];
    bzeroflag = Int32(param[6]);
    switchflag = Int32(param[7]);
    quadraticflag = Int32(param[8]);
    chemflag = Int32(param[9]);
    bnormflag = Int32(param[10]);
    wselfallflag = Int32(param[11]);   
    
    sna = initializesna(ntypes, twojmax, chemflag)

    facto = factable()    
    ccall(Libdl.dlsym(snaplib, :cpuInitSna), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint), 
            sna.rootpqarray, sna.cglist, facto, sna.idx_max, sna.idxz, sna.idxz_block, 
            sna.idxb, sna.idxb_block, sna.idxu_block, sna.idxcg_block, sna.twojmax)

    sna.radelem[2:end] = elemradius
    sna.wjelem[2:end] = elemweight    
    sna.map[2:end] = Int32.(Array(1:ntypes) .- 1)

    if snapcoeff === nothing       
        sna.coeffelem = zeros(Float64, sna.ncoeffall) 
    elseif length(snapcoeff) == sna.ntypes*sna.ncoeffall
        sna.coeffelem = snapcoeff
    else
        display([length(snapcoeff) sna.ntypes*sna.ncoeffall])
        error("number of snap coefficients is incompatible with the snap paramaters")
    end

    rcutmax = 0.0
    for ielem = 1:ntypes
        rcutmax = max(2.0*elemradius[ielem]*rcutfac,rcutmax)
    end

    for i = 1:ntypes
        for j = 1:ntypes
            cut = (elemradius[i] + elemradius[j])*rcutfac
            sna.rcutsq[j + (i-1)*(ntypes)] = cut*cut
        end
    end

    wself = 1.0;      
    www = wself*wself*wself;
    for j = 0:twojmax
        if (bnormflag == 1)
            sna.bzero[j+1] = www;
        else
            sna.bzero[j+1] = www*(j+1);
        end
    end

    sna.nperdim = sna.ncoeffall-1;
    sna.bnormflag = bnormflag;
    sna.chemflag = chemflag;    
    sna.quadraticflag = quadraticflag;
    sna.switchflag = switchflag;
    sna.bzeroflag = bzeroflag;
    sna.wselfallflag = wselfallflag;
    sna.wself = wself;
    sna.rmin0 = rmin0;
    sna.rfac0 = rfac0;
    sna.rcutfac = rcutfac;
    sna.rcutmax = rcutmax;        

    return sna
end

function snappotential(x, t, a, b, c, pbc, param, elemradius, elemweight, snapcoeff)

    sna = initsna(param, elemradius, elemweight, snapcoeff)
    eatom, fatom, vatom = snappotential(x, t, a, b, c, pbc, sna)
    return eatom, fatom, vatom 
end

function snappotential(x, t, a, b, c, pbc, sna)

    dim, N = size(x)
    
    #idxcg_max = sna.idxcg_max;
    idxu_max = sna.idxu_max;
    idxb_max = sna.idxb_max;
    idxz_max = sna.idxz_max;    
    twojmax = sna.twojmax;
    #ncoeff = sna.ncoeff;
    ncoeffall = sna.ncoeffall;
    ntypes = sna.ntypes;
    nelem = sna.nelements;    
    # ndoubles = sna.ndoubles;   
    # ntriples = sna.ntriples;   
    # nperdim = sna.nperdim;
    bnormflag = sna.bnormflag;
    chemflag = sna.chemflag;    
    quadraticflag = sna.quadraticflag;
    switchflag = sna.switchflag;    
    bzeroflag = sna.bzeroflag;
    wselfallflag = sna.wselfallflag;        

    map = sna.map;
    idxz = sna.idxz;
    #idxz_block = sna.idxz_block;
    idxb = sna.idxb;
    idxb_block = sna.idxb_block;
    idxu_block = sna.idxu_block;
    idxcg_block = sna.idxcg_block;   
    
    rcutmax = sna.rcutmax
    wself = sna.wself;
    rmin0 = sna.rmin0;
    rfac0 = sna.rfac0;
    rcutfac = sna.rcutfac;
    rcutmax = sna.rcutmax;        
    bzero = sna.bzero;
    rootpqarray = sna.rootpqarray;
    cglist = sna.cglist;
    rcutsq = sna.rcutsq;    
    radelem = sna.radelem;
    wjelem = sna.wjelem; 
    coeffelem = sna.coeffelem;           
    
    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:N));   
    atomtype = Int32.(t[:])
    alist = Int32.(alist)
    neighlist = Int32.(neighlist)
    neighnum = Int32.(neighnum)

    #pairlist, pairnum = neighpairlist(y, ilist, neighlist, neighnum, rcutmax*rcutmax);
    pairlist, pairnum = neighpairlist(y, ilist, alist, atomtype, neighlist, neighnum, rcutsq, ntypes)    
    rij, ai, aj, ti, tj = neighpairs(y, pairlist, pairnum, atomtype, ilist, alist);                
    ijnum = pairnum[end]
    
    # offset 1 for C++ code
    ilist = ilist .- Int32(1)
    alist = alist .- Int32(1)
    ai = ai .- Int32(1)
    aj = aj .- Int32(1)
    N = Int32(N)

    ulisttot_r = zeros(N*idxu_max*nelem)
    ulisttot_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeUi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
        map, ai, ti, tj, twojmax, idxu_max, N, ijnum, switchflag, chemflag)    

    ccall(Libdl.dlsym(snaplib, :cpuAddWself2Ui), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, wself, idxu_block, atomtype, map, wselfallflag, chemflag, 
        idxu_max, nelem, twojmax, N)    

    eatom = zeros(N)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeEi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        eatom, ulisttot_r, ulisttot_i, cglist, bzero, coeffelem, ilist, map, atomtype, 
        idxb, idxcg_block, idxu_block, twojmax, idxb_max, idxu_max, nelem, 
        ncoeffall, bnormflag, bzeroflag, wselfallflag, chemflag, N)    
        
    ylist_r = zeros(N*idxu_max*nelem)
    ylist_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeYi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, coeffelem, map, atomtype, idxz, idxb_block, 
        idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, nelem, ncoeffall, bnormflag, chemflag, N)    

    fatom = zeros(3,N)       
    vatom = zeros(6,N)            
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeFi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, 
        Cint, Cint, Cint, Cint), 
        fatom, vatom, ylist_r, ylist_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
        map, ai, ai, aj, ti, tj, twojmax, idxu_max, N, N, ijnum, switchflag, chemflag)    
    
    return eatom, fatom, vatom 
end

function snappotential(x, t, alist, neighlist, neighnum, sna::SnaStruct)

    dim = size(x,1);
    N = length(t)
    
    #idxcg_max = sna.idxcg_max;
    idxu_max = sna.idxu_max;
    idxb_max = sna.idxb_max;
    idxz_max = sna.idxz_max;    
    twojmax = sna.twojmax;
    #ncoeff = sna.ncoeff;
    ncoeffall = sna.ncoeffall;
    ntypes = sna.ntypes;
    nelem = sna.nelements;    
    # ndoubles = sna.ndoubles;   
    # ntriples = sna.ntriples;   
    # nperdim = sna.nperdim;
    bnormflag = sna.bnormflag;
    chemflag = sna.chemflag;    
    quadraticflag = sna.quadraticflag;
    switchflag = sna.switchflag;    
    bzeroflag = sna.bzeroflag;
    wselfallflag = sna.wselfallflag;        

    map = sna.map;
    idxz = sna.idxz;
    #idxz_block = sna.idxz_block;
    idxb = sna.idxb;
    idxb_block = sna.idxb_block;
    idxu_block = sna.idxu_block;
    idxcg_block = sna.idxcg_block;   
    
    rcutmax = sna.rcutmax
    wself = sna.wself;
    rmin0 = sna.rmin0;
    rfac0 = sna.rfac0;
    rcutfac = sna.rcutfac;
    rcutmax = sna.rcutmax;        
    bzero = sna.bzero;
    rootpqarray = sna.rootpqarray;
    cglist = sna.cglist;
    rcutsq = sna.rcutsq;    
    radelem = sna.radelem;
    wjelem = sna.wjelem; 
    coeffelem = sna.coeffelem;           
    
    #y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:N));   
    atomtype = Int32.(t[:])
    alist = Int32.(alist)
    neighlist = Int32.(neighlist)
    neighnum = Int32.(neighnum)

    pairlist, pairnum = neighpairlist(x, ilist, alist, atomtype, neighlist, neighnum, rcutsq, ntypes)    
    rij, ai, aj, ti, tj = neighpairs(x, pairlist, pairnum, atomtype, ilist, alist);                
    ijnum = pairnum[end]
    
    # offset 1 for C++ code
    ilist = ilist .- Int32(1)
    alist = alist .- Int32(1)
    ai = ai .- Int32(1)
    aj = aj .- Int32(1)

    ulisttot_r = zeros(N*idxu_max*nelem)
    ulisttot_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeUi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
        map, ai, ti, tj, twojmax, idxu_max, N, ijnum, switchflag, chemflag)    

    ccall(Libdl.dlsym(snaplib, :cpuAddWself2Ui), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, wself, idxu_block, atomtype, map, wselfallflag, chemflag, 
        idxu_max, nelem, twojmax, N)    

    eatom = zeros(N)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeEi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        eatom, ulisttot_r, ulisttot_i, cglist, bzero, coeffelem, ilist, map, atomtype, 
        idxb, idxcg_block, idxu_block, twojmax, idxb_max, idxu_max, nelem, 
        ncoeffall, bnormflag, bzeroflag, wselfallflag, chemflag, N)    
    
    ylist_r = zeros(N*idxu_max*nelem)
    ylist_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeYi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, coeffelem, map, atomtype, idxz, idxb_block, 
        idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, nelem, ncoeffall, bnormflag, chemflag, N)    

    fatom = zeros(3,N)       
    vatom = zeros(6,N)            
    ccall(Libdl.dlsym(snaplib, :cpuSnapComputeFi), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, 
        Cint, Cint, Cint, Cint), 
        fatom, vatom, ylist_r, ylist_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
        map, ai, ai, aj, ti, tj, twojmax, idxu_max, N, N, ijnum, switchflag, chemflag)    
        
    return eatom, fatom, vatom 
end

function snapdescriptors(x, t, a, b, c, pbc, param, elemradius, elemweight)
    sna = initsna(param, elemradius, elemweight)    
    return snapdescriptors(x, t, a, b, c, pbc, sna)
end

function snapdescriptors(x, t, a, b, c, pbc, snap::SNAPparams)
    sna = initsnap(snap)
    return  snapdescriptors(x, t, a, b, c, pbc, sna)
end

function snapdescriptors(x, t, a, b, c, pbc, sna::SnaStruct)

    dim, N = size(x)
    
    #idxcg_max = sna.idxcg_max;
    idxu_max = sna.idxu_max;
    idxb_max = sna.idxb_max;
    idxz_max = sna.idxz_max;    
    twojmax = sna.twojmax;
    ncoeff = sna.ncoeff;
    #ncoeffall = sna.ncoeffall;
    ntypes = sna.ntypes;
    nelem = sna.nelements;    
    ndoubles = sna.ndoubles;   
    ntriples = sna.ntriples;   
    #nperdim = sna.nperdim;
    bnormflag = sna.bnormflag;
    chemflag = sna.chemflag;    
    quadraticflag = sna.quadraticflag;
    switchflag = sna.switchflag;    
    bzeroflag = sna.bzeroflag;
    wselfallflag = sna.wselfallflag;        

    map = sna.map;
    idxz = sna.idxz;
    idxz_block = sna.idxz_block;
    idxb = sna.idxb;
    #idxb_block = sna.idxb_block;
    idxu_block = sna.idxu_block;
    idxcg_block = sna.idxcg_block;   
    
    rcutmax = sna.rcutmax
    wself = sna.wself;
    rmin0 = sna.rmin0;
    rfac0 = sna.rfac0;
    rcutfac = sna.rcutfac;
    rcutmax = sna.rcutmax;        
    bzero = sna.bzero;
    rootpqarray = sna.rootpqarray;
    cglist = sna.cglist;
    rcutsq = sna.rcutsq;    
    radelem = sna.radelem;
    wjelem = sna.wjelem; 
    #coeffelem = sna.coeffelem;                   

    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:N));   
    atomtype = Int32.(t[:])
    alist = Int32.(alist)
    neighlist = Int32.(neighlist)
    neighnum = Int32.(neighnum)
    
    pairlist, pairnum = neighpairlist(y, ilist, alist, atomtype, neighlist, neighnum, rcutsq, ntypes)    
    rij, ai, aj, ti, tj = neighpairs(y, pairlist, pairnum, atomtype, ilist, alist);                
    ijnum = pairnum[end]

    if ijnum>0
    # offset 1 for C++ code
    ilist = ilist .- Int32(1)
    alist = alist .- Int32(1)
    ai = ai .- Int32(1)
    aj = aj .- Int32(1)

    n = idxu_max*ijnum
    ulist_r = zeros(n)
    ulist_i = zeros(n)
    dulist_r = zeros(3*n)
    dulist_i = zeros(3*n)    
    ccall(Libdl.dlsym(snaplib, :cpuComputeUij), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, 
        Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
        ulist_r, ulist_i, dulist_r, dulist_i, rootpqarray, rij, wjelem, radelem, rmin0, 
        rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, ijnum, switchflag)    

    ulisttot_r = zeros(N*idxu_max*nelem)
    ulisttot_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuZeroUarraytot2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},  Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, wself, idxu_block, atomtype, map, ai, wselfallflag, chemflag, 
        idxu_max, nelem, twojmax, N)    

    ccall(Libdl.dlsym(snaplib, :cpuAddUarraytot), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, ulist_r, ulist_i, map, ai, tj, idxu_max, N, ijnum, chemflag)    

    zlist_r = zeros(idxz_max*ndoubles*N)
    zlist_i = zeros(idxz_max*ndoubles*N)
    ccall(Libdl.dlsym(snaplib, :cpuComputeZi2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, idxcg_block, twojmax, 
        idxu_max, idxz_max, nelem, bnormflag, N)    

    blist = zeros(idxb_max*ntriples*N)
    ccall(Libdl.dlsym(snaplib, :cpuComputeBi2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, ilist, atomtype, 
        map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
        nelem, bzeroflag, wselfallflag, chemflag, N)    
    
    dblist = zeros(idxb_max*ntriples*dim*ijnum)          
    ccall(Libdl.dlsym(snaplib, :cpuComputeDbidrj), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, idxz_block, map, ai, tj, 
        twojmax, idxb_max, idxu_max, idxz_max, nelem, bnormflag, chemflag, N, ijnum)    
    
    bi = zeros(ntypes*ncoeff)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrum), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},  
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint), 
        bi, blist, ilist, atomtype, N, ncoeff, ncoeff, ntypes, quadraticflag)    

    bd = zeros(ntypes*ncoeff*3*N)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrumDeriv), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},  
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        bd, blist, dblist, ai, ai, aj, ti, N, ijnum, ncoeff, ncoeff, ntypes, quadraticflag)    
    
    bv = zeros(6*ntypes*ncoeff)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrumVirial), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},   
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        bv, blist, dblist, rij, ai, ti, N, ijnum, ncoeff, ncoeff, ntypes, quadraticflag)  

    bd = reshape(bd, (3,N,ntypes*ncoeff))
    bv = reshape(bv, (6,ntypes*ncoeff))
    bv = bv[[1; 2; 3; 6; 5; 4],:]

    ai = ai .+ Int32(1)
    aj = aj .+ Int32(1)

    blist = reshape(blist, (N, idxb_max*ntriples)) 
    dblist = reshape(dblist, (ijnum, dim, idxb_max*ntriples)) 

    end

    if ijnum==0
    bi = zeros(ntypes*ncoeff)        
    bd = zeros(3,N,ntypes*ncoeff)
    bv = zeros(6,ntypes*ncoeff)    
    blist = zeros(N, idxb_max*ntriples)
    dblist = zeros(ijnum, dim, idxb_max*ntriples)    
    end

    return bi, bd, bv, blist, dblist, ai, aj, pairnum             
end


function snaponebodyefatom(x, t, a, b, c, pbc, descriptors)

    dim, N = size(x)
    M = length(descriptors.species)    
    eatom = zeros(N,M)
    for m = 1:M
        ind = findall(t[:] .== m);
        eatom[ind,m] .= 1.0
    end
    fatom = zeros(dim*N,M)

    return eatom, fatom 
end

function snaptwobodyefatom(x, t, a, b, c, pbc, descriptors)
    return nothing
end

function snapthreebodyefatom(x, t, a, b, c, pbc, descriptors)
    return nothing
end

function snapfourbodyefatom(x, t, a, b, c, pbc, sna)

    dim, N = size(x)
    
    #idxcg_max = sna.idxcg_max;
    idxu_max = sna.idxu_max;
    idxb_max = sna.idxb_max;
    idxz_max = sna.idxz_max;    
    twojmax = sna.twojmax;
    ncoeff = sna.ncoeff;
    #ncoeffall = sna.ncoeffall;
    ntypes = sna.ntypes;
    nelem = sna.nelements;    
    ndoubles = sna.ndoubles;   
    ntriples = sna.ntriples;   
    #nperdim = sna.nperdim;
    bnormflag = sna.bnormflag;
    chemflag = sna.chemflag;    
    quadraticflag = sna.quadraticflag;
    switchflag = sna.switchflag;    
    bzeroflag = sna.bzeroflag;
    wselfallflag = sna.wselfallflag;        

    map = sna.map;
    idxz = sna.idxz;
    idxz_block = sna.idxz_block;
    idxb = sna.idxb;
    #idxb_block = sna.idxb_block;
    idxu_block = sna.idxu_block;
    idxcg_block = sna.idxcg_block;   
    
    rcutmax = sna.rcutmax
    wself = sna.wself;
    rmin0 = sna.rmin0;
    rfac0 = sna.rfac0;
    rcutfac = sna.rcutfac;
    rcutmax = sna.rcutmax;        
    bzero = sna.bzero;
    rootpqarray = sna.rootpqarray;
    cglist = sna.cglist;
    rcutsq = sna.rcutsq;    
    radelem = sna.radelem;
    wjelem = sna.wjelem; 
    #coeffelem = sna.coeffelem;                   

    y, alist, neighlist, neighnum = fullneighborlist(x, a, b, c, pbc, rcutmax);
    ilist = Int32.(Array(1:N));   
    atomtype = Int32.(t[:])
    alist = Int32.(alist)
    neighlist = Int32.(neighlist)
    neighnum = Int32.(neighnum)
    
    pairlist, pairnum = neighpairlist(y, ilist, alist, atomtype, neighlist, neighnum, rcutsq, ntypes)    
    rij, ai, aj, ti, tj = neighpairs(y, pairlist, pairnum, atomtype, ilist, alist);                
    ijnum = pairnum[end]

    if ijnum>0
    # offset 1 for C++ code
    ilist = ilist .- Int32(1)
    alist = alist .- Int32(1)
    ai = ai .- Int32(1)
    aj = aj .- Int32(1)

    n = idxu_max*ijnum
    ulist_r = zeros(n)
    ulist_i = zeros(n)
    dulist_r = zeros(3*n)
    dulist_i = zeros(3*n)    
    ccall(Libdl.dlsym(snaplib, :cpuComputeUij), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, 
        Cdouble, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
        ulist_r, ulist_i, dulist_r, dulist_i, rootpqarray, rij, wjelem, radelem, rmin0, 
        rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, ijnum, switchflag)    

    ulisttot_r = zeros(N*idxu_max*nelem)
    ulisttot_i = zeros(N*idxu_max*nelem)
    ccall(Libdl.dlsym(snaplib, :cpuZeroUarraytot2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, 
        Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},  Cint, Cint, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, wself, idxu_block, atomtype, map, ai, wselfallflag, chemflag, 
        idxu_max, nelem, twojmax, N)    

    ccall(Libdl.dlsym(snaplib, :cpuAddUarraytot), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint), 
        ulisttot_r, ulisttot_i, ulist_r, ulist_i, map, ai, tj, idxu_max, N, ijnum, chemflag)    

    zlist_r = zeros(idxz_max*ndoubles*N)
    zlist_i = zeros(idxz_max*ndoubles*N)
    ccall(Libdl.dlsym(snaplib, :cpuComputeZi2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, idxz, idxu_block, idxcg_block, twojmax, 
        idxu_max, idxz_max, nelem, bnormflag, N)    

    blist = zeros(idxb_max*ntriples*N)
    ccall(Libdl.dlsym(snaplib, :cpuComputeBi2), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, bzero, ilist, atomtype, 
        map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
        nelem, bzeroflag, wselfallflag, chemflag, N)    
    
    dblist = zeros(idxb_max*ntriples*dim*ijnum)          
    ccall(Libdl.dlsym(snaplib, :cpuComputeDbidrj), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, 
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, 
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint, Cint), 
        dblist, zlist_r, zlist_i, dulist_r, dulist_i, idxb, idxu_block, idxz_block, map, ai, tj, 
        twojmax, idxb_max, idxu_max, idxz_max, nelem, bnormflag, chemflag, N, ijnum)    
    
    bi = zeros(ntypes*ncoeff)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrum), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},  
        Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint), 
        bi, blist, ilist, atomtype, N, ncoeff, ncoeff, ntypes, quadraticflag)    

    bd = zeros(ntypes*ncoeff*3*N)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrumDeriv), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},  
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        bd, blist, dblist, ai, ai, aj, ti, N, ijnum, ncoeff, ncoeff, ntypes, quadraticflag)    
    
    bv = zeros(6*ntypes*ncoeff)
    ccall(Libdl.dlsym(snaplib, :cpuSnapTallyBispectrumVirial), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},   
        Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint, Cint, Cint, Cint), 
        bv, blist, dblist, rij, ai, ti, N, ijnum, ncoeff, ncoeff, ntypes, quadraticflag)  

    bd = reshape(bd, (3*N,ntypes*ncoeff))
    bv = reshape(bv, (6,ntypes*ncoeff))
    bv = bv[[1; 2; 3; 6; 5; 4],:]

    ai = ai .+ Int32(1)
    aj = aj .+ Int32(1)

    blist = reshape(blist, (N, idxb_max*ntriples)) 
    dblist = reshape(dblist, (ijnum, dim, idxb_max*ntriples)) 

    end

    if ijnum==0
    bi = zeros(ntypes*ncoeff)        
    bd = zeros(3*N,ntypes*ncoeff)
    bv = zeros(6,ntypes*ncoeff)    
    blist = zeros(N, idxb_max*ntriples)
    dblist = zeros(ijnum, dim, idxb_max*ntriples)    
    end

    #return bi, bd, bv, blist, dblist, ai, aj, pairnum             

    eatom = bi
    fatom = bd 
    return eatom, fatom, blist 
end

function snapefatom(x, t, a, b, c, pbc, descriptors)

    globd = 0.0;
    fatom = 0.0;
    eatom1 = 0.0; 
    fatom1 = 0.0; 
    globd1 = 0.0; 
    for n = 1:length(descriptors)
        if (descriptors[n].name == "SNAP") & (descriptors[n].nbody==1) 
            eatom1, fatom1 = snaponebodyefatom(x, t, a[:], b[:], c[:], pbc[:], descriptors[n])                    
            globd1 = sum(eatom1, dims=1)            
        elseif (descriptors[n].name == "SNAP") & (descriptors[n].nbody==4) 
            sna = initsnap(descriptors[n])
            globd1, fatom1, eatom1 = snapfourbodyefatom(x, t, a[:], b[:], c[:], pbc[:], sna)                
            globd1 = reshape(globd1, (1, length(globd1)))         
            # display(x)
            # display(eatom1)
            # display(globd1)
            # display(sum(eatom1, dims=1))            
            # display(fatom1)
            # display(sna.bzeroflag)
            # display(sna.bnormflag)
            # error("here")
        end
        if (n==1) 
            globd = 1.0*globd1
            fatom = 1.0*fatom1 
        else                
            globd = cat(globd, globd1, dims=2)                
            fatom = cat(fatom, fatom1, dims=2)        
        end
    end            
    
    return globd, fatom    
end
    
function snaplinearsystem(config, descriptors, normalizeenergy)
    
    # cumalative sum of numbers of atoms 
    nconfigs = length(config.natom)
    natom = [0; cumsum(config.natom[:])];
    
    matA = 0.0
    vecb = 0.0
    for i = 1:nconfigs
        ci = i
        normconst = 1.0
        if normalizeenergy==1
            normconst = 1.0/config.natom[ci]
        end    
        x, t, a, b, c, pbc, e, f = getconfig(config, natom, ci)
        globd, fatom = snapefatom(x, t, a, b, c, pbc, descriptors)
        
        we2 = (config.we[1]*config.we[1])*(normconst*normconst)
        wf2 = (config.wf[1]*config.wf[1])
        matA = matA .+ (we2*(globd'*globd) + wf2*(fatom'*fatom))    
        vecb = vecb .+ ((we2*e)*(globd') + wf2*(fatom'*f[:]))    
    end
    
    return matA, vecb
    
end
    
function snapleastsq(data, descriptors, Doptions)
    
    display("Perform linear regression fitting ...")

    for i = length(descriptors):-1:1    
        if descriptors[i] === nothing
            deleteat!(descriptors, i)            
        end
    end   
    for i = length(data):-1:1    
        if data[i] === nothing
            deleteat!(data, i)            
        end
    end   

    pbc = Doptions.pbc
    normalizeenergy = Doptions.normalizeenergy
    a = Doptions.a
    b = Doptions.b
    c = Doptions.c

    n = length(data)
    matA = 0.0
    vecb = 0.0
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end
        A1, b1 = snaplinearsystem(config, descriptors, normalizeenergy)        
        matA = matA .+ A1
        vecb = vecb .+ b1         
    end    

    for i = 1:size(matA,1)
        matA[i,i] = matA[i,i]*(1.0 + 1e-12);
    end
    
    coeff = matA\vecb 
    return coeff, matA, vecb 
end
    
function snaperrors(config, descriptors, coeff, normalizeenergy)
    
    # cumalative sum of numbers of atoms 
    nconfigs = length(config.natom)
    natom = [0; cumsum(config.natom[:])];
    
    DFTenergy = zeros(nconfigs)
    energy = zeros(nconfigs)
    eerr = zeros(nconfigs)
    fmae = zeros(nconfigs)
    frmse = zeros(nconfigs)
    szbf = zeros(nconfigs)
    for i = 1:nconfigs
        ci = i
        normconst = 1.0
        if normalizeenergy==1
            normconst = 1.0/config.natom[ci]
        end    
        x, t, a, b, c, pbc, e, f = getconfig(config, natom, ci)
        globd, fatom = snapefatom(x, t, a, b, c, pbc, descriptors)
        e1 = globd*coeff;
        f1 = fatom*coeff;

        DFTenergy[i] = e*normconst      
        energy[i] = e1[1]*normconst        
        eerr[i] = abs(e1[1] - e)*normconst         

        N = length(f1)
        szbf[i] = N 
        res = (f1 - f[:])    
        fmae[i] = sum(abs.(res))/N     
        ssr = sum(res.^2)
        mse = ssr /N;
        frmse[i] = sqrt(mse) 
    end
    
    emae = sum(abs.(eerr))/nconfigs
    ssr = sum(eerr.^2)
    mse = ssr/nconfigs
    ermse = sqrt(mse) 

    nforces = sum(szbf)
    fmaeave = sum(fmae.*szbf)/nforces
    frmseave =  sqrt(sum((frmse.^2).*szbf)/nforces)    

    return eerr, emae, ermse, fmae, frmse, fmaeave, frmseave, nconfigs, nforces, energy, DFTenergy
end
        
function snaperroranalysis(data, descriptors, Doptions, coeff)
    
    display("Calculate errors ...")

    for i = length(descriptors):-1:1    
        if descriptors[i] === nothing
            deleteat!(descriptors, i)            
        end
    end   
    for i = length(data):-1:1    
        if data[i] === nothing
            deleteat!(data, i)            
        end
    end   

    pbc = Doptions.pbc
    normalizeenergy = Doptions.normalizeenergy
    a = Doptions.a
    b = Doptions.b
    c = Doptions.c

    n = length(data)
    emae = -ones(Float64,n)
    fmae = -ones(Float64,n)
    ermse = -ones(Float64,n)
    frmse = -ones(Float64,n)
    szbe = zeros(Int64, n)
    szbf = zeros(Int64, n)
    DFTenergies = Array{Any}(nothing, n)
    energies = Array{Any}(nothing, n)
    eerr = Array{Any}(nothing, n)
    ferr = Array{Any}(nothing, n)
    ferm = Array{Any}(nothing, n)
    for i = 1:n
        if typeof(data[i]) == Preprocessing.DataStruct
            config, ~ = Preprocessing.readconfigdata(data[i], pbc, a, b, c)                  
        else
            config = data[i]
        end

        e1, e2, e3, e4, e5, e6, e7, nconfigs, nforces, energy, DFTenergy = 
                snaperrors(config, descriptors, coeff, normalizeenergy)        

        DFTenergies[i] = DFTenergy                 
        energies[i] = energy                 
        eerr[i] = e1
        ferr[i] = e4 
        ferm[i] = e5 

        szbe[i] = nconfigs
        szbf[i] = nforces
        emae[i] = e2
        ermse[i] = e3 
        fmae[i] = e6
        frmse[i] = e7 
    end    
    
    energyerrors = zeros((n+1),2)    
    forceerrors = zeros((n+1),2)    
    stresserrors = zeros((n+1),2)    

    energyerrors[1,1] = sum(emae.*szbe)/sum(szbe)  # MAE
    energyerrors[1,2] = sqrt(sum((ermse.^2).*szbe)/sum(szbe)) # RMSE 
    energyerrors[2:end,:] = [emae ermse]        

    forceerrors[1,1] = sum(fmae.*szbf)/sum(szbf)
    forceerrors[1,2] =  sqrt(sum((frmse.^2).*szbf)/sum(szbf))    
    forceerrors[2:end,:] = [fmae frmse]   

    return energyerrors, forceerrors, stresserrors, eerr, ferr, ferm, energies, DFTenergies  
end
    
    