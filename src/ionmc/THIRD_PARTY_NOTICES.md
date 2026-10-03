# Third-party notices

## Geant4

Parts of this repository derive from Geant4 source code, release 11.4.2:

- the tabulated Ashley-Ritchie-Brandt function used for the Barkas correction in
  `src/ionmc/physics/stopping.py` and the structure of the Barkas, Bloch and Mott terms follow
  `G4EmCorrections.cc`;
- the ICRU 90 parser in `src/ionmc/data/icru90.py` expects the array layout of
  `G4ICRU90StoppingData.cc` (the file itself is downloaded by users, not redistributed here;
  small excerpts appear as test fixtures);
- material compositions, densities and mean excitation energies in `src/ionmc/materials.py`
  are taken from `G4NistMaterialBuilder.cc`.

This product includes software developed by Members of the Geant4 Collaboration
( http://cern.ch/geant4 ).

The Geant4 Software License follows verbatim, as published in the file `LICENSE` of the
Geant4 v11.4.2 source tree (https://github.com/Geant4/geant4/blob/v11.4.2/LICENSE). The
web page https://geant4.web.cern.ch/download/license only lists the copyright holders
and refers to this license.

```text
Geant4 Software License
Version 1.0, 28 June 2006

Copyright (c) Copyright Holders of the Geant4 Collaboration, 1994-2006.
See http://cern.ch/geant4/license for details on the copyright holders.  All
rights not expressly granted under this license are reserved.

This software includes voluntary contributions made to Geant4.
See http://cern.ch/geant4 for more information on Geant4.

Installation, use, reproduction, display, modification and redistribution of
this software, with or without modification, in source and binary forms, are
permitted on a non-exclusive basis. Any exercise of rights by you under this
license is subject to the following conditions:

1. Redistributions of this software,  in whole or in part,  with  or without
   modification, must reproduce the above copyright notice and these license
   conditions  in  this  software,  the  user  documentation  and  any other
   materials provided with the redistributed software.

2. The user documentation,if any,included with a redistribution,must include
   the following notice:"This product includes software developed by Members
   of the Geant4 Collaboration ( http://cern.ch/geant4 )."
   If that  is  where  third-party  acknowledgments  normally  appear,  this
   acknowledgment  must  be  reproduced  in  the  modified  version  of this
   software itself.

3. The names "Geant4" and "The Geant4 toolkit" may not be used to endorse or
   promote software,or products derived therefrom, except with prior written
   permission  by  license@geant4.org.  If this software is redistributed in
   modified form,  the name  and  reference of  the modified version must be
   clearly distinguishable from that of this software.

4. You are under no obligation to provide anyone with  any  modifications of
   this software that you may develop,including but not limited to bug fixes,
   patches,  upgrades  or  other enhancements or derivatives of the features,
   functionality or performance of this software. However, if you publish or
   distribute your modifications without  contemporaneously  requiring users
   to enter into a separate written license agreement,  then  you are deemed
   to have granted all  Members  and all  Copyright Holders  of  the  Geant4
   Collaboration  a license to your  modifications,  including modifications
   protected by any patent owned by you,under the conditions of this license.

5. You may not include  this  software in  whole or in part in any patent or
   patent  application  in  respect  of  any  modification  of this software
   developed by you.


6. DISCLAIMER

 THIS SOFTWARE IS PROVIDED BY THE MEMBERS AND COPYRIGHT HOLDERS OF THE GEANT4
 COLLABORATION AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 INCLUDING,  BUT NOT LIMITED TO,  IMPLIED WARRANTIES OF  MERCHANTABILITY,  OF
 SATISFACTORY QUALITY,  AND FITNESS  FOR  A  PARTICULAR PURPOSE  OR  USE  ARE
 DISCLAIMED. THE MEMBERS OF THE GEANT4 COLLABORATION AND CONTRIBUTORS MAKE NO
 REPRESENTATION THAT THE SOFTWARE AND MODIFICATIONS THEREOF,WILL NOT INFRINGE
 ANY PATENT, COPYRIGHT, TRADE SECRET OR OTHER PROPRIETARY RIGHT.

7. LIMITATION OF LIABILITY

 THE  MEMBERS  AND  COPYRIGHT   HOLDERS  OF  THE   GEANT4  COLLABORATION  AND
 CONTRIBUTORS SHALL HAVE NO LIABILITY FOR DIRECT,INDIRECT,SPECIAL, INCIDENTAL,
 CONSEQUENTIAL,  EXEMPLARY,  OR PUNITIVE DAMAGES  OF  ANY CHARACTER INCLUDING,
 WITHOUT LIMITATION, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES, LOSS OF USE,
 DATA OR PROFITS,  OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND ON ANY THEORY
 OF CONTRACT,  WARRANTY, TORT  (INCLUDING NEGLIGENCE),  PRODUCT LIABILITY  OR
 OTHERWISE,  ARISING IN  ANY WAY  OUT OF THE USE  OF  THIS SOFTWARE,  EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.

8. This license  shall terminate with immediate effect and without notice if
   you  fail  to  comply with any  of  the terms of this license,  or if you
   institute litigation against any Member or Copyright Holder of the Geant4
   Collaboration with regard to this software.
```

## NIST Standard Reference Database 124 (ESTAR, PSTAR, ASTAR)

Copyright protection on this compilation of data has been secured by the Secretary of the
U.S. Department of Commerce on behalf of the United States.

M.J. Berger, J.S. Coursey, M.A. Zucker, J. Chang (2005), ESTAR, PSTAR, and ASTAR: Computer
Programs for Calculating Stopping-Power and Range Tables for Electrons, Protons, and Helium
Ions (version 2.0.1), National Institute of Standards and Technology, Gaithersburg, MD,
https://doi.org/10.18434/T4NC7P. Terms: https://www.nist.gov/open/license.

The NIST tables are downloaded by each user through NIST's public web service into a
local cache and are not shipped with this package. The current source tree contains no
NIST table values or values from which they can be recovered; see decision 0038 in the
repository for the project's legal assessment and for an inventory of earlier repository
history that quoted individual values.


## ICRU Report 90

International Commission on Radiation Units and Measurements, Report 90: Key Data for
Ionizing-Radiation Dosimetry: Measurement Standards and Applications, Journal of the ICRU
14(1) (2014, published 2016). The report is copyrighted and is not redistributed here. The
proton and alpha stopping-power values used for evaluation are those embedded in the Geant4
source file named above.
