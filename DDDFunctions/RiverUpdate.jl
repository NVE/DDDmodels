#    Function River_update
#
#------------------------------------------------------------------------
#     Description:  Shifts to on timestep ahead on the QRivx from previous timestep and updates water in the river network
#
#     Author: Thomas Skaugen
#     Revised: 16.12.2019
#--------------------------------------------------------------------------

function RiverUpdate!(vectorlengde, QRivx, qsimutx)
# qsimutx: 1 dim array float
# QRivx: 1 dim array float
# QRD : float    
# vectorlengde: scalar integer

    #vectorlengde = 10 
    #QRivx = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
    #qsimutx = [3.2, 3.1, 3.0, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3,2.2]
    
  if vectorlengde > 1
    @views QRivx[1:(vectorlengde-1)] .= QRivx[2:vectorlengde]
    QRivx[vectorlengde] = 0.0
    QRivx .+= qsimutx
  else
     QRivx .= qsimutx
  end  
  QRD = QRivx[1] # QRD is Runoff for current timestep

return QRD
end            
