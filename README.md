## Relative Energy Calculation  

This is a simple MD engine. Different dihedrals of C-C-C-C are generated first. Then, non-bonded and bonded interactions for each snapshot are calculated. It should be mentioned that non-bonded interactions are treated in a way that is 1-4 and with a scaling factor of 0.5, scaled Coulomb and LJ potentials. For the coulomb term, the relative dielectric constant is considered as 1.  
In the case of bonded interactions, harmonic potentials are used for bonds and angles, while Ryckaert-Bellemans function is used for proper dihedrals.


    ![Image](https://github.com/user-attachments/assets/4706641c-4f4b-4ef1-8a2e-4ab0b062bd35)
