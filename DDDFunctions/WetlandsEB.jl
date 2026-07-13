#    Function Wetlands(isoilx,sncovx,smx,outx,ceax,tempx,lp,ep)
#
#C------------------------------------------------------------------------
#C
#      New moisture routine 
#C
#C     Beskrivelse:  Rutinen beregner markvannsmagasin. 
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#
#     Parametere:
#     ea       Actual evapotransp   
#     toSoil    Rain/snowmelt
#     smbog       Soil moisture content, mm # note that smbog is the entire reservoir for bogs, so delta_smbog is dicharge to the RN
#     D       Saturation deficit (D>Z) mm
#     M       Total subsurface water reservoir mm
#     OUTBOG  Outflow form soil moisture reservoir mm/d
#--------------------------------------------------------------------------


function WetlandsEB(misoil,eatemp,middelsca,smbog,M,ET)

#-------------------------------------------------------------------------
#  Actual evapotranspiration
#-------------------------------------------------------------------------

       if (eatemp > -10.0)              # Function of temp, areal  meanm, see Sælthun 1996, p. 9                          
         eabog = min(ET, ET*((smbog+misoil)/M))        
       else          
         eabog = 0.0
       end
      
#------------------------------------------------------------------------
#     Updating sm caused by evapotranspiration
#------------------------------------------------------------------------
      smbog -= eabog
      
      if (smbog < 0) # prevents breaching the water balance, must adjust ea      
        eabog += smbog #sm is negative
        smbog = 0
      end
      
#------------------------------------------------------------------------
#     Updating smbog caused by evapotranspiration
#------------------------------------------------------------------------
      Bograt = (smbog+misoil)
    
      if(Bograt > M)      
        outbog = (Bograt-M)    #Excess water released to runoff
        smbog = M      
      else      
        outbog = 0
        smbog += misoil
      end 
                                            
return outbog, smbog, eabog
end


