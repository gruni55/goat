2022-05-30
<p>
  New branch for elastic scattering added. 
  Phasejump detection for LightSrcPlane (tubedray) now working. 
  </p> <br>
2022-09-20
<p>
  All classes of GOAT are now encapsulated into namespace GOAT. <br>
  Sub-namespaces: <br>
  <ul>
    <li><b>maths (GOAT::maths)</b> is used for all mathematical classes, functions and constants (everything inside directory "maths") </li>
    <li><b>raytracing (GOAT::raytracing)</b> is used for all classes, functions and constants related to the raytracing part (everything inside directory "raytracing) </li>
    </ul>
</p>
</p> <br>
2022-12-17 Classes for inelastic scattering added. Now it is also possible to calculate and store the electric field inside an object. Also new light
sources with arbitrary ray distribution is provided for gaussian and for plane beam. In the moment, necessary normalisation to the total power is not done,
i.e. a given power is not concerned yet. 
<p>

