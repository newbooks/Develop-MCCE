#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"

using namespace std;

struct ccoi //coi contains circle of intersection data for pairs when objects are involved
{
    SGrid <delphi_real> xyz; //vector applied in the center of the sphere and pointing to the center of the coi
    delphi_real rad;         //coi radius;
    delphi_integer is;       //number of the intersections between an atom and an object
};

void CDelphiSpace::sas()
{
    delphi_integer nver=520,nedge=1040;

    SGrid <delphi_real> ver[nver+1];
    delphi_integer edgv[2+1][nedge+1],edg[nedge+1],oti[nver+1],st[nedge+1];

    //a third field in pls takes into account the correspondence nprobj-npr

    delphi_integer jprec,nprtobj,nprobj;
    delphi_integer ii,jj,dim,dim1;
    string strtmp,strtmp1;

    SGrid <delphi_real> xyzm,x123, dx123, tij123;
    SGrid <delphi_real> rmv[3+1],cf123, dy123;
    SGrid <delphi_integer> ix123, ic123;
    vector < SGrid <delphi_integer> > pls; // pls is changed to be local variable
    vector < SGrid <delphi_integer> > plstemp;

    ccoi coi_1;
    vector <ccoi> coi;
    vector <ccoi> coitemp;
    vector < SGrid <delphi_real> > expostemp;

    delphi_real rad2;
    delphi_integer nacc,nacct,i,j,k,ie,ia2,ie2,ie1,ilvl,ia1,iv,iv1,iv2;
    delphi_integer ip,liml,limu,ne,nlvl,npr,nprp,nprx,nst,nprt,nvo;
    delphi_integer nxa,nvi, nv;
    delphi_real ctf,ctf2,cst,dctf,d2,del,dx1,dx2,dx3,ds2,pre,rad, radj;
    delphi_real rij,rdn,rv1,rv2,rvmg,dctf2,sm1,sm2,snt,tta,vmg;
    delphi_real cbln,csp,dmg,tm;

    coi_1.xyz= sgrid_temp_real;
    coi_1.rad=0.;
    coi_1.is=0;

    fill_n(ver,nver+1,sgrid_temp_real);

#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif

    for (i=0; i<3; i++)
    {
        for (j=0; j<nedge+1; j++)
        {
            edgv[i][j]=0;
        }
    }

    nacc=extot;
    nacct=0;
    nprt=0;
    nprtobj=0;
    nlvl=5;
    nvi=12;

    radpmax=max(fRadPrb[1],fRadPrb[2]);

    cbln=2.*(fRMax+radpmax); //cube length

    cubedata(1.0,cbln);

    dim=(lcb+1)*(mcb+1)*(ncb+1); // dim: how many cubes

    cbn1_v.assign(dim+1,0);
    cbn2_v.assign(dim+1,0);

    dim1=27;

    if ((iNObject-numbmol)>0) dim1=max(dim,delphi_integer(27));

    cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);

    cube();

    tta=2.*fPi/nvi; //generate a template of vertices....

#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif

    for(i=1; i<=nvi; i++)
    {
        rdn=(i-1)*tta;

        ver[i].nX=cos(rdn);
        ver[i].nY=sin(rdn);
        ver[i].nZ=0.;

        j=i+1;
        if(i==nvi)j=1;
        edgv[1][i]=i;
        edgv[2][i]=j;
    }

    nv=nvi;
    ne=nvi;
    ie1=1;
    ie2=0;

    for(ilvl=1; ilvl<=nlvl; ilvl++)
    {
        ie1=ie2+1;
        ie2=ne;
        for(ie=ie1; ie<=ie2; ie++)
        {
            iv1=edgv[1][ie];
            iv2=edgv[2][ie];

            //2011-05-26 Using operations on coord and int_coord type
            //variables defined in module operators_on_coordinates
            xyzm=ver[iv1]+ver[iv2];
            vmg=sqrt(optDot(xyzm,xyzm));

            nv=nv+1;
            ver[nv]=xyzm/vmg;
            ne=ne+1;
            edg[ie]=ne;
            edgv[1][ne]=iv1;
            edgv[2][ne]=nv;
            ne=ne+1;
            edgv[1][ne]=nv;
            edgv[2][ne]=iv2;
        }
    }

    ne=ie2;

#ifdef PARALLEL_OMP
#pragma omp for schedule(auto)
#endif

    for (i=ie1; i<=ne; i++)
    {
        edg[i]=-1;
    }

    #ifdef VERBOSE
    cout <<" # of vertices             :" << setw(10) << nv << endl
            <<" # of edges                :" << setw(10) << ne << endl;
    #endif

    for (i=0; i<iNatom+1; i++)
    {
        ast[i]=1;
    }

    nacc=0;
    npr=0;
    nprobj=0;
    nprp=0;

    //it finds pairs.......
    for(i=1; i<=iNatom; i++)
    {
        rad=r0[i];
        rad2=r02[i];
        if(sDelPhiPDB[i].radius==0.) continue;
        x123=xn1[i];
        ix123=optCast <delphi_integer,delphi_real> ((x123-xyzo)*cbai);

        liml=cbn1[ix123.nX][ix123.nY][ix123.nZ];
        limu=cbn2[ix123.nX][ix123.nY][ix123.nZ];

        //2011-06-17 Resizing of arrays keeping old values intact
        if ((npr+limu-liml+1)>nprt)
        {
            nprt=nprt+5000;

            if(pls.size() > 0)
                pls.resize(nprt+1);
            else
                pls.assign(nprt+1, sgrid_temp_int);
        }

        if ((nprobj+limu-liml+1)>nprtobj)
        {
            nprtobj=nprtobj+1000;

            if(coi.size() >0)
                coi.resize(nprobj+1);
            else
                coi.assign(nprobj+1,coi_1);
        }

        jprec=0;
        for(jj=liml; jj<=limu; jj++)
        {
            j=cbal[jj];

            if (j<=iNatom)
            {
                radj=r0[j];

                if (sDelPhiPDB[j].radius>0.&&j>i)
                {
                    ctf=rad+radj;
                    ctf2=ctf*ctf;
                    dctf=abs(rad-radj);
                    dctf2=dctf*dctf;
                    dx123=xn1[j]-x123;
                    d2=optDot(dx123,dx123);
                    del=ctf2-d2;

                    if (del>0.01&&d2>dctf2)
                    {
                        npr=npr+1;
                        pls[npr].nX=i;
                        pls[npr].nY=j;
                        pls[npr].nZ=0;
                    }
                }
            }
            else //objects are abandoned:
            {

            }

            jprec=j;

        }

        if (npr==nprp) ast[i]=0;

        nprp=npr;
    }

    if(cbn1_v.size()>0) vector <delphi_integer> ().swap(cbn1_v);
    if(cbn2_v.size()>0) vector <delphi_integer> ().swap(cbn2_v);
    if(cbal.size()>0) vector <delphi_integer> ().swap(cbal);

    if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
    if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);

    #ifdef VERBOSE
    cout <<" # of pairs                :" << setw(10) << npr << endl;
    #endif

    //2011-05-26 Temporarily removed time calculations
    cbln=fRMax+radpmax;

    cubedata(2.0,cbln);

    dim1=27;
    if ((iNObject-numbmol)>0) dim1=max(dim, delphi_integer(27));

    cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);
    cube();

    nprx=0;
    for(ip=1; ip<=npr; ip++)
    {
        i=pls[ip].nX;
        j=pls[ip].nY;

        if (j<=iNatom)
        {
            dx123=xn1[j]-xn1[i];
            d2=optDot(dx123,dx123);
            dmg=sqrt(d2);
            pre=1.+(r02[i]-r02[j])/d2;
            tij123=xn1[i]+((0.5*pre)*dx123);
            rij=0.5 * sqrt((r0[i]+r0[j])*(r0[i]+r0[j])-d2) * sqrt(d2-(r0[i]-r0[j])*(r0[i]-r0[j]))/dmg;
        }
        else // never go to this else statement
        {
            cout << "### Warning: if j < iNatom else:" << endl;
            nprobj=pls[ip].nZ;
            rij=coi[nprobj].rad;

            //pay attention, here dx has a different meaning from previous one
            dx123=coi[nprobj].xyz;
            d2=optDot(dx123,dx123);
            tij123=xn1[i]+dx123;
            dmg=sqrt(d2);
        }

        dx1=dx123.nX;
        dx2=dx123.nY;
        dx3=dx123.nZ;
        rvmg=sqrt(dx1*dx1+dx2*dx2);

        if (rvmg>1.0e-8)
        {
            rv1=-dx2/rvmg;
            rv2=dx1/rvmg;
            cst=dx3/dmg;

            snt=sqrt(1.-cst*cst); //snt=rvmg/dmg !doesn't lead to any improved performance

            csp=1.0-cst;
            tm=csp*rv1;
            sm1=snt*rv1;
            sm2=snt*rv2;

            rmv[1].nX=tm*rv1+cst; //rmv[1]= {tm*rv1+cst,tm*rv2,sm2};
            rmv[1].nY=tm*rv2;
            rmv[1].nZ=sm2;

            rmv[2].nX=tm*rv2; //rmv[2]= {tm*rv2,csp*rv2*rv2+cst,-sm1};
            rmv[2].nY=csp*rv2*rv2+cst;
            rmv[2].nZ=-sm1;

            rmv[3].nX=-sm2; //rmv[3]= {-sm2,sm1,cst};
            rmv[3].nY=sm1;
            rmv[3].nZ=cst;
        }
        else
        {
            rmv[1].nX=1.0; //rmv[1]= {1.,0.,0.};
            rmv[1].nY=0.0;
            rmv[1].nZ=0.0;

            rmv[2].nX=0.0; //rmv[2]= {0.,1.,0.};
            rmv[2].nY=1.0;
            rmv[2].nZ=0.0;

            rmv[3].nX=0.0; //rmv[3]= {0.,0.,1.};
            rmv[3].nY=0.0;
            rmv[3].nZ=1.0;
        }

        nvo=0;

        /*
         * assign memory to expos if needed
         */
        //2011-06-17 Re-sizing array keeping old value
        if ((nacc+nv)>nacct)
        {
            nacct=nacct+1000;

            if (expos.size()>0)
                expos.resize(nacct+1);
            else
                expos.assign(nacct+1,sgrid_temp_real);
        }

        iv=1;

        D10:
        while(iv<=nvi)
        {
            /*
             * +rm(7)*ver[3][iv] has been removed because it is always zero
             */
            cf123.nX=rmv[1].nX*ver[iv].nX+rmv[1].nY*ver[iv].nY;
            cf123.nY=rmv[2].nX*ver[iv].nX+rmv[2].nY*ver[iv].nY;
            cf123.nZ=rmv[3].nX*ver[iv].nX+rmv[3].nY*ver[iv].nY;
            cf123=tij123+(cf123*rij);

            ic123=optCast <delphi_integer,delphi_real> ((cf123-xyzo)*cbai);

            liml=cbn1[ic123.nX][ic123.nY][ic123.nZ];
            limu=cbn2[ic123.nX][ic123.nY][ic123.nZ];


            ii=liml;

            D05:
            while(ii<=limu)
            {
                k=cbal[ii];
                if (k>iNatom)
                {
                    oti[iv]=k;
                    ii++;
                    goto D05;
                }

                dy123=xn1[k]-cf123;
                ds2=optDot(dy123,dy123);

                if (ds2<rs2[k])
                {
                    oti[iv]=k;
                    iv++;
                    goto D10;
                }

                ii++;
            }

            nvo=nvo+1;
            nacc=nacc+1;
            expos[nacc]=cf123;

            oti[iv]=0;

            iv++;
        }

        nst=0;
        if (nlvl>0)
        {
            for(ie=nvi; ie>=1; ie--)
            {
                ia1=oti[edgv[1][ie]];
                ia2=oti[edgv[2][ie]];

                if(ia1>0&&ia1==ia2) continue;
                nst=nst+1;
                st[nst]=ie;
            }
        }

        if (nst>0)
        {
            D030:
            while(true)
            {
                ie=st[nst];
                nst=nst-1;
                ia1=oti[edgv[1][ie]];
                ia2=oti[edgv[2][ie]];

                if ((ia1>iNatom)||(ia2>iNatom))
                {
                    if(nst>0) goto D030;
                    break;
                }

                iv=ie+nvi;

                cf123.nX=optDot(rmv[1],ver[iv]);
                cf123.nY=optDot(rmv[2],ver[iv]);
                cf123.nZ=optDot(rmv[3],ver[iv]);
                cf123=tij123+(cf123*rij);

                if (ia1!=0)
                {
                    dy123=xn1[ia1]-cf123;
                    ds2=optDot(dy123,dy123);

                    if (ds2<rs2[ia1])
                    {
                        oti[iv]=ia1;

                        if (edg[ie]>0)
                        {
                            nst=nst+1;
                            st[nst]=edg[ie]+1;
                        }

                        if(nst>0) goto D030;

                        break;
                    }
                }

                if (ia2!=0)
                {
                    dy123=xn1[ia2]-cf123;
                    ds2=optDot(dy123,dy123);

                    if (ds2<rs2[ia2])
                    {
                        oti[iv]=ia2;

                        if (edg[ie]>0)
                        {
                            nst=nst+1;
                            st[nst]=edg[ie];
                        }

                        if(nst>0) goto D030;

                        break;
                    }
                }

                ic123=optCast <delphi_integer,delphi_real> ((cf123-xyzo)*cbai);
                liml=cbn1[ic123.nX][ic123.nY][ic123.nZ];
                limu=cbn2[ic123.nX][ic123.nY][ic123.nZ];

                ii=liml;

                D055:
                while(ii<=limu)
                {
                    k=cbal[ii];

                    if (k>iNatom)
                    {
                        oti[iv]=k;
                        ii++;
                        goto D055;
                    }// if

                    dy123=xn1[k]-cf123;
                    ds2=optDot(dy123,dy123);

                    if (ds2<rs2[k])
                    {
                        oti[iv]=k;

                        if (edg[ie]>0)
                        {
                            nst=nst+1;
                            st[nst]=edg[ie]+1;
                            nst=nst+1;
                            st[nst]=edg[ie];
                        }

                        if(nst>0) goto D030;
                        goto END030;
                    }
                    ii++;
                }

                nvo=nvo+1;
                nacc=nacc+1;
                expos[nacc]=cf123;

                oti[iv]=0;
                if (edg[ie]>0)
                {
                    if ( edg[edg[ie]+1] > 0 || ia2>0 )
                    {
                        nst=nst+1;
                        st[nst]=edg[ie]+1;
                    }

                    if (edg[edg[ie]]>0||ia1>0)
                    {
                        nst=nst+1;
                        st[nst]=edg[ie];
                    }
                }

                if(nst<=0) goto END030;
            }
            END030:

            cout << "";

        }

        if (nvo>0) //considering pairs also where one 'partner' is an object
        {
            nprx=nprx+1;
            ast[i]=0;
            if (j<=iNatom) ast[j]=0;
        }
    }

    if(pls.size()>0)    vector < SGrid <delphi_integer> >().swap(pls);
    if(coi.size()>0)    vector <ccoi> ().swap(coi);
    if(cbn1_v.size()>0) vector <delphi_integer> ().swap(cbn1_v);
    if(cbn2_v.size()>0) vector <delphi_integer> ().swap(cbn2_v);
    if(cbal.size()>0)   vector <delphi_integer> ().swap(cbal);

    if(cbn1 != NULL) free_pt3d<delphi_integer>(cbn1,lcb+1,mcb+1,ncb+1);
    if(cbn2 != NULL) free_pt3d<delphi_integer>(cbn2,lcb+1,mcb+1,ncb+1);


    nxa=0;
    for(i=1; i<=iNatom; i++)
        if(ast[i]==0) nxa=nxa+1;

    #ifdef VERBOSE
    cout <<"# pairs analyzed (atom-atom and atom-object)= " << npr << endl;
    cout <<"# exposed pairs (atom-atom and atom-object)= " << nprx << endl;
    cout <<"no. arc points = " << nacc << endl;
    cout <<"no. surface atoms = " << nxa << " nbur = " << iNatom-nxa << endl;
    #endif

    extot=nacc;

}// void sas;
