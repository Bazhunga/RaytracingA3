>> Original
Render
   Shoots out rays
   shadeRay
      traverseScene
         does a bunch of intersections
      checks if intersection occurs
         if yes computesShading
            Iterates thorugh lightsource
            shades(ray)
               does the specular stuff


>> Modified

Render
   Shoots out rays
   shadeRay
      traverseScene
         does a bunch of intersections
      checks if intersection occurs
         if yes, shoot more rays (will need a counter to keep track of the iteration #)
         computesShading
            Iterates thorugh lightsource
            shades(ray)
               does the specular stuff

