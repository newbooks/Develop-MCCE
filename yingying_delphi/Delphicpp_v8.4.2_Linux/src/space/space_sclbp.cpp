#include "../interface/environment.h"
#include "../delphi/delphi_constants.h"
#include "../misc/misc_grid.h"
#include "../space/space.h"

using namespace std;


/**
 * program to reposition the boundary grid points on the molecular surface
 * (S. Sridharan May 1994)
 */
void CDelphiSpace::sclbp()
{
    delphi_integer nbra[SPACE_NBRA_SIZE];
    bool out;
    bool outcb[5][5][5]; //outcb index should be modified

    delphi_integer iaprec,dim1,kind;
    delphi_integer dim,prevmed,med;
    delphi_integer imezzo[7],iente[7],ix,iy,iz;
    SGrid <delphi_integer > ixyz;
    delphi_integer iac1,iac2;
    bool lga,lgd,iflag,precedenza,vicinanza,flag;
    string strtmp;
    SGrid <delphi_real> vnor,s123;
    delphi_real radpmax,dst;
    SGrid <delphi_real> x1xyz,xg, u123,xxyyzz;
    delphi_real temp;
    SGrid <delphi_real> xq,dxyz,dixyz,dx123, dr;
    delphi_real hgs,ds2min1,ds2min2;
    SGrid <delphi_integer > it,jxyz;

    delphi_integer iac,ia,i,j,k,iacl,iii,ii,jjx,jjy,jjz,jx,jy,jz;
    delphi_integer jzi,jyi,jxi,liml,limu,kk,nnbr,ncbp;
    delphi_real x1,cbln,del,dis2,dis,dist,dmn,dcr,ctf,dmn1,dmn2;
    delphi_real cba,ds2,dsr,rmn,rdist,dmx;

    vector <bool> internal;

    internal.assign(iNObject+1,false);
    //epsdim=iNatom+iNObject+2; replaced in space module semi global
    iflag=false;
    iac1=0;
    iac2=0;

    iall=0;

    /*
     * hgs= half grid spacing
     */
    hgs=1./(2.*fScale);

    for (i=0; i<=4; i++) //outcb=true; initiallized above.
        for (j=0; j<=4; j++)
            for (k=0; k<=4; k++)
                outcb[i][j][k]=true;

    for (i=1; i<=3; i++) //outcb[-1:1][-1:1][-1:1]=false;
        for (j=1; j<=3; j++)
            for (k=1; k<=3; k++)
                outcb[i][j][k]=false;

    x1=1.0/fScale; // conversion from grid to delphi_real coordinates(can also use routine gtoc)

    //x1xyz=cOldMid-(0.5*x1*optCast<delphi_real,delphi_integer>(global_iGrid+1)); //for parallel
    if (ibctyp == 3)
        x1xyz = cOldMid - double( 0.5*x1*(iGrid.nX + 1) );
    else
        x1xyz = myStartCoor - 0.5 *x1;

    radpmax=max(fRadPrb[1],fRadPrb[2]);

    if (extot==0&&radpmax>0.0&&(iNObject>1||iNatom>1))
    {
        //find extrema
        //here one should consider the global system (Walter)
        #ifdef VERBOSE
        cout <<"Scaling routine in action//" << endl;
        #endif

        //2011-05-19 Using operations on coord and int_coord type variables defined in module operators_on_coordinates
        cMin.nX=6000.;
        cMin.nY=6000.;
        cMin.nZ=6000.;

        cMax.nX=-6000.;
        cMax.nY=-6000.;
        cMax.nZ=-6000.;

        for(ii=0; ii<iNObject; ii++)
        {
            cMin=optMin(cMin,sLimObject[ii].nMin);
            cMax=optMax(cMax,sLimObject[ii].nMax);
        }

        sas();
    }

    del=radpmax;
    del=max(del,1./(2.*fScale));
    cbln=fRMax+del;

    cubedata(2.0,cbln);

    dim1=27;
    if ((iNObject-numbmol)>0) dim1=max(dim, delphi_integer(27));
    cbal.assign(dim1*(iNatom+iNObject-numbmol)+1,0);

    cube();

    ncbp=0;

    #ifdef VERBOSE
    cout << "iBoundNum: " << iBoundNum << endl;
    #endif

    for(i=1; i<=iBoundNum; i++)
    {
        //+per trattare molecole con diversa epsilon++01/2002+
        if (iBoundNum!=iBoundNumsurf&&numbmol>1)
        {
            #ifdef VERBOSE
            cout << "###### sclbp iente iepsmp: " << endl;
            #endif

            ixyz=optCast <delphi_integer,delphi_real> (scspos[i]);

            ix=ixyz.nX;
            iy=ixyz.nY;
            iz=ixyz.nZ;

            iflag=false;

            iente[1]=iepsmp[ix][iy][iz].nX%epsdim;
            iente[2]=iepsmp[ix][iy][iz].nY%epsdim;
            iente[3]=iepsmp[ix][iy][iz].nZ%epsdim;
            iente[4]=iepsmp[ix-1][iy][iz].nX%epsdim;
            iente[5]=iepsmp[ix][iy-1][iz].nY%epsdim;
            iente[6]=iepsmp[ix][iy][iz-1].nX%epsdim;

            imezzo[1]=iepsmp[ix][iy][iz].nX/epsdim;
            imezzo[2]=iepsmp[ix][iy][iz].nY/epsdim;
            imezzo[3]=iepsmp[ix][iy][iz].nZ/epsdim;
            imezzo[4]=iepsmp[ix-1][iy][iz].nX/epsdim;
            imezzo[5]=iepsmp[ix][iy-1][iz].nY/epsdim;
            imezzo[6]=iepsmp[ix][iy][iz-1].nZ/epsdim;

            //guardo se ho due molecole con diversa epsilon nel punto,interno
            if(imezzo[1]!=imezzo[6]&&imezzo[1]*imezzo[6]!=0) iflag=(iente[1]<=iNatom+1&&iente[6]<=iNatom+1);

            //iflag sar?vero se il bgp ?interno ma non con un oggetto
            for(ii=2; ii<=6; ii++)
                if(imezzo[ii]!=imezzo[ii-1] && imezzo[ii]*imezzo[ii-1]!=0)
                    iflag=iflag||(iente[ii]<=iNatom+1 && iente[ii-1]<=iNatom+1);
        }

        xg=scspos[i]*x1+x1xyz;

        /*
         * find the closest surface atom to the gridpoint
         */
        it=optCast <delphi_integer,delphi_real> ((xg-xyzo)*cbai);

        dmn=100.;
        ds2min1=1000.;
        ds2min2=1000.;
        prevmed=0;
        iac=0;
        nnbr=0;

        lmncb.nX=lcb;
        lmncb.nY=mcb;
        lmncb.nZ=ncb;

        if (optORLT(it, delphi_integer(0))||optORGT(it,lmncb))
        {
            // if the bgp is outside the cube, probably it is due to some object
            for(ii=1; ii<=iNObject; ii++)
            {
                strtmp=dataobject_v[(ii-1)*2];
                kind = atoi(strtmp.substr(15,3).c_str());

                if (strtmp.substr(0,4)!="is a" && kind!=2)
                {
                    if ( optANDLE(xg,(sLimObject[ii].nMax+x1) )&& optANDGT(xg,(sLimObject[ii].nMin-x1)) )
                    {
                        nnbr=nnbr+1;
                        if(nnbr >= SPACE_NBRA_SIZE)
                            cout << "space_sclbp>> index beyond size of nbra: nnbr= "<< nnbr << endl;
			nbra[nnbr]=ii+iNatom;
                        liml=1;
                        limu=0;
                    }
                }
            }

            #ifdef VERBOSE
            if(liml!=1||limu!=0) cout <<"bgp close to nothing" << endl;
            #endif
        }
        else
        {
            liml=cbn1[it.nX][it.nY][it.nZ];
            limu=cbn2[it.nX][it.nY][it.nZ];
        }

        iaprec=0;

        for(kk=liml; kk<=limu; kk++)
        {
            ia=cbal[kk];

            if (ia<=iNatom)
            {
                // added iflag to still save in atsurf value
                if (iflag)
                {
                    dx123=xg-xn1[ia];
                    dis2=optDot(dx123,dx123)-sDelPhiPDB[ia].radius*sDelPhiPDB[ia].radius;

                    // dis2, and ds2min are distance**2 from center - radius**2 so they can be <0
                    if (dis2<ds2min1)
                    {
                        iac2=iac1;
                        ds2min2=ds2min1;
                        iac1=ia;
                        ds2min1=dis2;
                    }
                    else if (dis2<=ds2min2)
                    {
                        iac2=ia;
                        ds2min2=dis2;
                    }
                }
                else
                {
                    if (ast[ia]==0)
                    {
                        nnbr=nnbr+1;
                        if(nnbr >= SPACE_NBRA_SIZE)
                            cout << "space_sclbp>> index beyond size of nbra: nnbr= "<< nnbr << endl;
                        nbra[nnbr]=ia;
                    }
                }
            }
            else
            {
                if (ia!=iaprec)
                {
                    iaprec=ia;
                    nnbr=nnbr+1;
                    if(nnbr >= SPACE_NBRA_SIZE)
                        cout << "space_sclbp>> index beyond size of nbra: nnbr= "<< nnbr << endl;
                    nbra[nnbr]=ia;
                }
            }
        }

        if (iflag)
        {
            atsurf[i]=iac1;

            if (iac1*iac2==0||iac1==iac2)
            {
                cout <<"Problems in Scaling multidielectric Boundary Grid Points" << endl;
                exit(0);
            }

            atndx[i]=-1;

            // looks like atndx is used to build Delunay surface, so excluding these bgps
            dx123=xn1[iac2]-xn1[iac1];
            temp=optDot(dx123,dx123);
            temp=0.5*(ds2min2-ds2min1)/temp;
            scspos[i]=xg+(temp*dx123);
            scsnor[i]=sgrid_temp_real;

            continue;
        }
        else
        {
            for(ii=1; ii<=nnbr; ii++)
            {
                if(ii >= SPACE_NBRA_SIZE)
                    cout << "space_sclbp>> index beyond size of nbra: ii= "<< ii << endl;
		ia=nbra[ii];
                med=iAtomMed[ia];
                lgd=(med!=prevmed);

                if (ia>iNatom)
                {
                    iii=ia-iNatom;
                    xq=xg;

                    // try to find closest VdW surface, better if it is buried
                    // internal is used for the object to which surface the bgp is closer
                    precedenza=ia>iac&&(iac>iNatom||iac==0);
                    vicinanza=abs(dist)<abs(dmn);
                    lga=(precedenza&&(vicinanza||dist<0.))||(vicinanza&&dmn>0.);

                    if ((dist<dmn&&!lgd)||(lga&&lgd))
                    {
                        dmn=dist;
                        iac=ia;
                        prevmed=med;
                        dr=dixyz*(dist-fRadPrb[1]);
                        vnor=dixyz;
                    }

                    internal[iii]=(dist<0.0);
                }
                else
                {
                    dx123=xg-xn1[ia];
                    dis=sqrt(optDot(dx123,dx123) )-sDelPhiPDB[ia].radius;
                    precedenza=ia>iac||iac>iNatom;
                    vicinanza=abs(dis)<abs(dmn);
                    lga=(precedenza&&(vicinanza||dis<0.))||(vicinanza&&dmn>0.);

                    if ((dis<dmn&&!lgd)||(lga&&lgd))
                    {
                        prevmed=med;
                        dmn=dis;
                        iac=ia;
                    }
                }
            }

            atsurf[i]=iac;
        }

        if (iac==0&&iac1==0)
        {
            cout <<"no close atom or object for boundary pointeger " << i << endl;
            exit(0);
        }

        // if iac is an object dr has alredy been calculated and HAS a DIFFERENT value!!!!!!
        if (iac<=iNatom) dr=xg-xn1[iac];

        dsr=sqrt( optDot(dr,dr));
        out=true;

        if (radpmax>0.0)
        {
            //u should have the same value as previous one
            if (iac<=iNatom)
                u123=xn1[iac]+(((r0[iac]*dr)/dsr));
            else
                u123=xg-dr;

            it=optCast <delphi_integer,delphi_real> ((u123-xyzo)*cbai);
            nnbr=0;

            liml=cbn1[it.nX][it.nY][it.nZ];
            limu=cbn2[it.nX][it.nY][it.nZ];

            for(kk=liml; kk<=limu; kk++)
            {
                ia=cbal[kk];
                if (ia<=iNatom)
                {
                    dx123=u123-xn1[ia];
                    ds2=optDot(dx123,dx123);
                    if(ds2<rs2[ia])out=false;
                }
                else
                {
                    if (ia!=iac&&(!internal[ia-iNatom]))
                    {
                        xq=u123; //I want to know if u is within the shell surrounding the object

                        if (dist>0.0&&dist<fRadPrb[1]-1.e-6) out=false;
                    }
                }
            }
        }

        if (out)
        {
            ncbp=ncbp+1;
            if (iac<=iNatom)
            {
                scspos[i]=xn1[iac]+(dr*(sDelPhiPDB[iac].radius/dsr));
                scsnor[i]=dr/dsr;
            }
            else
            {
                scspos[i]=(xg-(fRadPrb[1]*vnor))-dr;
                scsnor[i]=vnor;

            }

            atndx[i]=iac;
        }
        else
        {
            atndx[i]=0;
        }
    }

    //fScale the re-entrant points with respect to expos if fRadPrb = 0.0 we are done.
    if (radpmax>0.0)
    {
        iall=0;
        cba=1./grdi;

        for (i=1; i<=iBoundNum; i++)
        {
            //b+++mol. con diversa eps ++++01/02+++++++++++++++
            if(atndx[i]==-1) continue;

            if (atndx[i]==0)
            {
                s123=scspos[i]*x1+x1xyz;

                xxyyzz=(s123-mnxyz)*grdi;
                jxyz=optCast <delphi_integer,delphi_real> (xxyyzz);
                jx=jxyz.nX;
                jy=jxyz.nY;
                jz=jxyz.nZ;

                dxyz=xxyyzz-optCast<delphi_real,delphi_integer>(jxyz);

                dmn1=min(dxyz.nX,min(dxyz.nY,dxyz.nZ));
                dmx=max(dxyz.nX,max(dxyz.nY,dxyz.nZ));
                dmn2=1.0-dmx;
                dcr=min(dmn1,dmn2);
                ctf=cba*(1+dcr);

                ctf=ctf*ctf;
                iacl=0;
                rmn=100.;
                for(jjx=jx-1; jjx<=jx+1; jjx++)
                {
                    for(jjy=jy-1; jjy<=jy+1; jjy++)
                    {
                        for(jjz=jz-1; jjz<=jz+1; jjz++)
                        {
                            for(ii=iab1[jjx][jjy][jjz]; ii<=iab2[jjx][jjy][jjz]; ii++)
                            {
                                iac= icume[ii];
                                dist=optDot( (s123-expos[iac]),(s123-expos[iac]) );

                                if (dist<rmn)
                                {
                                    rmn=dist;
                                    iacl=iac;
                                }
                            }
                        }
                    }
                }

                if (!(iacl>0&&rmn<ctf))
                {
                    for(jxi=-2; jxi<=2; jxi++)
                    {
                        for(jyi=-2; jyi<=2; jyi++)
                        {
                            for(jzi=-2; jzi<=2; jzi++)
                            {
                                if (outcb[jxi+2][jyi+2][jzi+2]) // index of outcb is modified
                                {
                                    jjx=jx+jxi;
                                    if (jjx>=0&&jjx<=lcb1)
                                    {
                                        jjy=jy+jyi;
                                        if (jjy>=0&&jjy<=mcb1)
                                        {
                                            jjz=jz+jzi;
                                            if (jjz>=0&&jjz<=ncb1)
                                            {
                                                for(ii=iab1[jjx][jjy][jjz]; ii<=iab2[jjx][jjy][jjz]; ii++)
                                                {
                                                    iac= icume[ii];
                                                    dist= optDot( (s123-expos[iac]),(s123-expos[iac]) );

                                                    if (dist<rmn)
                                                    {
                                                        rmn=dist;
                                                        iacl=iac;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (iacl<=0)
                    {
                        iall=iall+1;
                        for(iac=1; iac<=extot; iac++)
                        {
                            dist=optDot((s123-expos[iac]),(s123-expos[iac]));
                            if (dist<rmn)
                            {
                                rmn=dist;
                                iacl=iac;
                            }
                        }
                    }
                }

                dxyz=s123-expos[iacl];
                rdist=sqrt(optDot(dxyz,dxyz) );

                if (rdist==0)
                {
                    dist=0.0;
                }
                else
                {
                    //if inside any object  fRadPrb[2]...
                    dst=0.;
                    flag=true;

                    for(ii=1; ii<=iNObject; ii++)
                    {
                        strtmp=dataobject_v[(ii-1)*2];
                        kind = atoi(strtmp.substr(15,3).c_str());

                        if (strtmp.substr(0,4)!="is a" && kind!=2)
                        {
                            xq=s123;

                            //assuming that if the VdW pointeger is half grid space into an object that means that this belongs to an atom buried in the object
                            if (dst<-hgs)
                            {
                                dist=fRadPrb[2]/rdist;
                                flag=false;
                                break;
                            }
                        }
                    }

                    if(flag) dist=fRadPrb[1]/rdist;
                }

                scspos[i]=expos[iacl]+(dxyz*dist);

                if (rdist>1.0e-8)
                {
                    scsnor[i]=(-dxyz)/rdist;
                }
                else
                {
                    #ifdef VERBOSE
                    cout <<"bdp close to arcp " << i << rdist << endl;
                    #endif
                }
            }
        }
    }

    if(internal.size()>0) vector <bool>().swap(internal);

}// void sclbp;

