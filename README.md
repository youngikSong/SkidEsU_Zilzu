# SkidEsU_Zilzu
 SkidEsU (Skidpad Estimation Utility) is Matlab code developed by KAIST Zilzu FSAE team, to predict vehicle's performance on skidpad and cornering in order to support engineers optimize their design and setup.

Skidpad : Fix the constraints of the car, with skidpad radius R. Code will generate MMM diagram at constant R

MMM : Basically same with Skidpad, but input is V instead of R. Calculates performance of the car at constant velocity

Development History - 

2024-04-25 : skidpad and mmm are first developed
2024-04-27 : Minor fix and plot improvement

Future Work -

SkidEsU is not complete! It currently doesn't calculate vehicle's roll steer effect, 3d geometry effect, roll center change, cog height changing by pitching, camber change due to longitudinal acceleration, and so on! I hope everybody to share their knowledge and ability to improve the code, so we can learn and help every engineers of FSAE to develop their team.

Roadmap - 

1. Camber change by longitudinal acceleration implementation
2. Non-constant roll center implementation
3. Roll steer implementation
4. 3d geometry implementation by entering coordinates of suspension points
5. COG height changing by pitch implementation
6. GUI implementation