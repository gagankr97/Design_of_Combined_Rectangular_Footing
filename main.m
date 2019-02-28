% Design of Combined Footing Footing
clc
clear
format short g

load input.mat
load value.mat
disp ("Design of Combined Rectangular Footing \n")

%Compute the Area of Footing 
Vertical_Load_on_Column_A = Load_on_Column_A;
printf("Vertical_Load_on_Column_A = %d KN \n", Vertical_Load_on_Column_A)

Vertical_Load_on_Column_B = Load_on_Column_B;
printf("Vertical_Load_on_Column_B = %d KN \n", Vertical_Load_on_Column_B)

Assume_the_Self_Weight_of_Footing = Self_Weight*100;
printf("Assume_the_Self_Weight_of_Footing = %d %% \n", Assume_the_Self_Weight_of_Footing)

Self_Weight_of_the_Footing = Self_Weight*(Vertical_Load_on_Column_A+Vertical_Load_on_Column_B);
printf("Self_Weight_of_the_Footing = %d KN \n", Self_Weight_of_the_Footing)

Total_Load = Vertical_Load_on_Column_A + Vertical_Load_on_Column_B + Self_Weight_of_the_Footing; 
printf("Total_Load = %d KN \n", Total_Load)
disp("\n")

Required_Area_of_Footing =  (Total_Load/Soil_Pressure);
printf("Required_Area_of_Footing = %d m^2 \n", Required_Area_of_Footing)
 
Provided_Area_of_Footing =round(Total_Load/Soil_Pressure);
printf("Provided_Area_of_Footing = %d m^2 \n", Provided_Area_of_Footing)

% Side of Rectangular Footing
Length_of_Longer_Side_of_Footing = round((sqrt(3*(Provided_Area_of_Footing))));
printf("Length_of_Longer_Side_of_Footing = %d m \n", Length_of_Longer_Side_of_Footing)
Length_of_Shorter_Side_of_Footing = (Provided_Area_of_Footing/Length_of_Longer_Side_of_Footing);
printf("Length_of_Shorter_Side_of_Footing = %d m \n", Length_of_Shorter_Side_of_Footing)
disp("\n")

Center_to_Center_spacing_of_Column = Center_to_Center_spacing_of_Column;
printf("Center_to_Center_spacing_of_Columns = %d m \n", Center_to_Center_spacing_of_Column)


xbar = round((Load_on_Column_B*Center_to_Center_spacing_of_Column)/(Load_on_Column_A+Load_on_Column_B));
Dist_btw_Side_of_the_footing_to_center_of_column_A = (((Length_of_Longer_Side_of_Footing)/2)-xbar);
printf("Dist_btw_Side_of_the_footing_to_center_of_column_A = %d m \n", Dist_btw_Side_of_the_footing_to_center_of_column_A)
Dist_btw_Side_of_the_footing_to_center_of_column_B = (Length_of_Longer_Side_of_Footing)-(Center_to_Center_spacing_of_Column+Dist_btw_Side_of_the_footing_to_center_of_column_A);
printf("Dist_btw_Side_of_the_footing_to_center_of_column_B = %d m \n", Dist_btw_Side_of_the_footing_to_center_of_column_B)
disp("\n")

Net_upward_pressure = (Load_on_Column_A + Load_on_Column_B)/(Length_of_Longer_Side_of_Footing * Length_of_Shorter_Side_of_Footing);
printf("Net_upward_pressure = %d KN/m^2 \n", Net_upward_pressure)
Pressure_per_meter_length = Net_upward_pressure * Length_of_Shorter_Side_of_Footing;
printf("Pressure_per_meter_length = %d KN/m^2 \n", Pressure_per_meter_length)
disp("\n")

% Shear Force
Shear_Force_just_left_of_A = Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A; 
printf("Shear_Force_just_left_of_A = %d KN \n", Shear_Force_just_left_of_A)
Shear_Force_just_right_of_A = Load_on_Column_A-Shear_Force_just_left_of_A;
printf("Shear_Force_just_right_of_A = %d KN \n", Shear_Force_just_right_of_A)
Shear_Force_just_right_of_B = Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_B;
printf("Shear_Force_just_right_of_B = %d KN \n", Shear_Force_just_right_of_B)
Shear_Force_just_left_of_B = Load_on_Column_B-Shear_Force_just_right_of_B;
printf("Shear_Force_just_left_of_B = %d KN \n", Shear_Force_just_left_of_B)
Shear_Force_is_zero_from_A_at_the_dist = Shear_Force_just_right_of_A/Pressure_per_meter_length;
printf("Shear_Force_is_zero_from_A_at_the_dist = %d m \n", Shear_Force_is_zero_from_A_at_the_dist)
disp("\n")

% Bending Moment
Maximum_Bending_Moment = round(((Load_on_Column_A * Shear_Force_is_zero_from_A_at_the_dist)-((Pressure_per_meter_length/2)*(Dist_btw_Side_of_the_footing_to_center_of_column_A + Shear_Force_is_zero_from_A_at_the_dist)^2)));
printf("Maximum_Bending_Moment = %d KNm \n", Maximum_Bending_Moment)
Sagging_Bending_Moment_at_A = ((Pressure_per_meter_length/2)*(Dist_btw_Side_of_the_footing_to_center_of_column_A)^2);
printf("Sagging_Bending_Moment_at_A = %d KNm \n", Sagging_Bending_Moment_at_A)
Sagging_Bending_Moment_at_B = ((Pressure_per_meter_length/2)*(Dist_btw_Side_of_the_footing_to_center_of_column_B)^2);
printf("Sagging_Bending_Moment_at_B = %d KNm \n", Sagging_Bending_Moment_at_B)
dist_bwt_side_of_footing_and_side_of_column_A = Dist_btw_Side_of_the_footing_to_center_of_column_A - (Width1_of_Column_A/2);
dist_bwt_side_of_footing_and_side_of_column_B = Dist_btw_Side_of_the_footing_to_center_of_column_B - (Width1_of_Column_B/2);

Sagging_Bending_Moment_at_outer_face_of_Column_A = round(((Pressure_per_meter_length)*(dist_bwt_side_of_footing_and_side_of_column_A)^2)/2);
printf("Sagging_Bending_Moment_at_outer_face_of_Column_A = %d KNm \n", Sagging_Bending_Moment_at_outer_face_of_Column_A)
Sagging_Bending_Moment_at_outer_face_of_Column_B = round(((Pressure_per_meter_length)*(dist_bwt_side_of_footing_and_side_of_column_B)^2)/2);
printf("Sagging_Bending_Moment_at_outer_face_of_Column_B = %d KNm \n", Sagging_Bending_Moment_at_outer_face_of_Column_B)

disp("\n")


x = ((-(Load_on_Column_A-Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A))+sqrt((Load_on_Column_A-Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A)*(Load_on_Column_A-Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A)...
-Pressure_per_meter_length*Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A*Dist_btw_Side_of_the_footing_to_center_of_column_A))/(-Pressure_per_meter_length);

x1 = ((-(Load_on_Column_A-Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A))-sqrt((Load_on_Column_A-Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A)*(Load_on_Column_A-Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A)...
-Pressure_per_meter_length*Pressure_per_meter_length*Dist_btw_Side_of_the_footing_to_center_of_column_A*Dist_btw_Side_of_the_footing_to_center_of_column_A))/(-Pressure_per_meter_length);

x2 = Center_to_Center_spacing_of_Column-x1;

Shear_Force_at_First_Point_of_Contraflexure = Shear_Force_just_right_of_A-(Pressure_per_meter_length*x);
printf("Shear_Force_at_First_Point_of_Contraflexure = %d KN \n", Shear_Force_at_First_Point_of_Contraflexure)
Shear_Force_at_Second_Point_of_Contraflexure = Shear_Force_just_left_of_B-(Pressure_per_meter_length*x2);
printf("Shear_Force_at_Second_Point_of_Contraflexure = %d KN \n", Shear_Force_at_Second_Point_of_Contraflexure)
disp("\n")

% Effective Depth from Bending Compression

xu_max_by_d = (700)/(1100+0.87*Fy);
Ru = (0.36*Fck*xu_max_by_d*(1-0.416*xu_max_by_d));


Effective_Depth = round((sqrt((1.5*Maximum_Bending_Moment*1000000)/(Ru*Length_of_Shorter_Side_of_Footing*1000)))/100)*100;
printf("Effective_Depth = %d mm \n", Effective_Depth)
%Check for punching shear
Width_of_critical_plane = (Width1_of_Column_A+Width1_of_Column_B);
Punching_Shear_Force = (1.5*(Load_on_Column_B-(Net_upward_pressure*Width_of_critical_plane*Width_of_critical_plane)));
printf("Punching_Shear_Force = %d KN \n", Punching_Shear_Force)
Actual_Shear_Stress = ((Punching_Shear_Force*1000)/(4*Width_of_critical_plane*1000*Effective_Depth));
printf("Actual_Shear_Stress = %d N/mm^2 \n", Actual_Shear_Stress)
Permissible_Shear_Stress = 0.25*(sqrt(Fck));
printf("Permissible_Shear_Stress = %d N/mm^2 \n", Permissible_Shear_Stress)
if (Actual_Shear_Stress>Permissible_Shear_Stress)
Required_Effective_Depth = ceil((Punching_Shear_Force*1000)/(4*Width_of_critical_plane*1000*Permissible_Shear_Stress));
printf("Required_Effective_Depth = %d mm \n", Required_Effective_Depth)
elseif(Actual_Shear_Stress<Permissible_Shear_Stress)
Required_Effective_Depth = Effective_Depth;
printf("Required_Effective_Depth = %d mm \n", Required_Effective_Depth)
endif

Overall_Depth = round((Required_Effective_Depth+Clear_Cover)/10)*10;
printf("Overall_Depth = %d mm \n", Overall_Depth)

Provided_Effective_Depth = Overall_Depth-Clear_Cover;
printf("Provided_Effective_Depth = %d mm \n", Provided_Effective_Depth)


disp("\n")

% Design for Bending Tension in Longitudinal direction
% Reinforcement for Hogging Bending Moment
disp("Reinforcement for Hogging Bending Moment")
Factored_Moment = 1.5*Maximum_Bending_Moment*1000000;
printf("Factored_Moment = %d Nmm \n", Factored_Moment)
Area_of_steel = 0.5*(Fck/Fy)*(1-sqrt(1-((4.6*Factored_Moment)/(Fck*Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth*Required_Effective_Depth))))*Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth;
printf("Area_of_steel = %d mm^2 \n", Area_of_steel)

Area_of_one_bar = (pi/4)*(dia*dia);
printf("Area_of_one_bar = %d mm^2 \n", Area_of_one_bar)
No_of_bars = ceil(Area_of_steel/Area_of_one_bar)

Provided_area_of_steel = No_of_bars*Area_of_one_bar;
printf("Provided_area_of_steel = %d mm^2 \n", Provided_area_of_steel)
Percentage_of_Steel = (Provided_area_of_steel*100)/(Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth);
printf("Percentage_of_Steel = %d %% \n", Percentage_of_Steel)

tbd = interp1 (table(:,1),table(:,2),Fck);
Ld = (Fy*dia)/(4*tbd*1.6);
printf("Development_Lenght = %d mm \n", Ld)
Ast = Provided_area_of_steel;
xu = (0.87*Fy*Ast)/(0.36*Fck*Length_of_Shorter_Side_of_Footing*1000);
M1 = 0.87*Fy*Ast*(Required_Effective_Depth-0.416*xu);

V1 = 1.5 * Shear_Force_at_First_Point_of_Contraflexure;
d = Required_Effective_Depth;
d1 = 12*dia;
if (d>d1)
  Lo = d;
elseif (d<d1)
  Lo = d1;
endif


if (((M1/V1)+Lo)>Ld)
  disp("Hence Code Requirements are Satisfied")
elseif 
  disp("Hemce Code Requirements arenot Satisfied")
endif

disp("\n")

% Reinforcement for Sagging Bending Moment
disp("Reinforcement for Sagging Bending Moment")
Factored_Moment_S1 = 1.5*Sagging_Bending_Moment_at_outer_face_of_Column_B*1000000;
printf("Factored_Moment = %d Nmm \n", Factored_Moment_S1)
Area_of_steel_2 = 0.5*(Fck/Fy)*(1-sqrt(1-((4.6*Factored_Moment_S1)/(Fck*Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth*Required_Effective_Depth))))*Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth;
printf("Area_of_steel = %d mm^2 \n", Area_of_steel_2)
Area_of_one_bar_2 = (pi/4)*(dia_12*dia_12);
printf("Area_of_one_bar = %d mm^2 \n", Area_of_one_bar_2);
No_of_bars_2 = ceil(Area_of_steel_2/Area_of_one_bar_2);
printf("No_of_bars = %d \n", No_of_bars_2)
Provided_area_of_steel_2 = No_of_bars_2*Area_of_one_bar_2;
printf("Provided_area_of_steel = %d mm^2 \n", Provided_area_of_steel_2)
Percentage_of_Steel_2 = (Provided_area_of_steel_2*100)/(Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth);
printf("Percentage_of_Steel = %d %% \n", Percentage_of_Steel_2)
tbd = interp1 (table(:,1),table(:,2),Fck);
Ld_2 = (Fy*dia_12)/(4*tbd*1.6);

printf("Development_Lenght = %d mm \n", Ld_2)
Ast_2 = Provided_area_of_steel_2;
xu_2 = (0.87*Fy*Ast_2)/(0.36*Fck*Length_of_Shorter_Side_of_Footing*1000);
M2 = 0.87*Fy*Ast_2*(Required_Effective_Depth-0.416*xu_2);

V2 = 1.5 * Shear_Force_at_First_Point_of_Contraflexure;
d = Required_Effective_Depth;
d2 = 12*dia_12;
if (d>d2)
  Lo = d;
elseif (d<d2)
  Lo = d2;
endif


if (((M2/V2)+Lo)>Ld_2)
  disp("Hence Code Requirements are Satisfied")
elseif 
  disp("Hemce Code Requirements arenot Satisfied")
endif

disp("\n")

% Reinforcement for Sagging Bending Moment
disp("Reinforcement for Sagging Bending Moment")
Factored_Moment_S2 = 1.5*Sagging_Bending_Moment_at_outer_face_of_Column_A*1000000;
printf("Factored_Moment = %d Nmm \n", Factored_Moment_S2)
Area_of_steel_3 = 0.5*(Fck/Fy)*(1-sqrt(1-((4.6*Factored_Moment_S2)/(Fck*Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth*Required_Effective_Depth))))*Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth;
printf("Area_of_steel = %d mm^2 \n", Area_of_steel_3)
Area_of_one_bar_3 = (pi/4)*(dia_12*dia_12);
printf("Area_of_steel = %d mm^2 \n", Area_of_one_bar_3)
if (Area_of_steel_3<((0.12*Overall_Depth*Length_of_Shorter_Side_of_Footing*1000)/100))
Provided_Area_of_steel_3=((0.12*Overall_Depth*Length_of_Shorter_Side_of_Footing*1000)/100);
printf("Provided_Area_of_steel = %d mm^2 \n", Provided_Area_of_steel_3)
endif

No_of_bars_3 = ceil(Provided_Area_of_steel_3/Area_of_one_bar_3);
printf("No_of_bars = %d \n", No_of_bars_3)
Percentage_of_Steel_2 = (Provided_Area_of_steel_3*100)/(Length_of_Shorter_Side_of_Footing*1000*Required_Effective_Depth);
printf("Percentage_of_Steel = %d %% \n", Percentage_of_Steel_2)


tbd = interp1 (table(:,1),table(:,2),Fck);
Ld_3 = (Fy*dia_12)/(4*tbd*1.6);

printf("Development_Lenght = %d mm \n", Ld_3)
Ast_3 = Provided_Area_of_steel_3;
xu_3 = (0.87*Fy*Ast_3)/(0.36*Fck*Length_of_Shorter_Side_of_Footing*1000);
M3 = 0.87*Fy*Ast_3*(Required_Effective_Depth-0.416*xu_3);

V3 = 1.5 * Shear_Force_at_First_Point_of_Contraflexure;
d = Required_Effective_Depth;
d3 = 12*dia_12;
if (d>d3)
  Lo = d;
elseif (d<d3)
  Lo = d3;
endif


if (((M3/V3)+Lo)>Ld_3)
  disp("Hence Code Requirements are Satisfied")
elseif 
  disp("Hemce Code Requirements arenot Satisfied")
endif

disp("\n")

% Check for One Way Shear
disp("Check for Shear")
disp("Check for One Way Shear Diagonal Tension from the centre of the Column B")
d = Provided_Effective_Depth/1000;
Cantilever_portion_B = d + (Width1_of_Column_B/2);
Shear_Force_1 = 1.5*(Shear_Force_just_right_of_B-(Pressure_per_meter_length*Cantilever_portion_B));
printf("Shear_Force = %d KN \n", Shear_Force_1)

tv1 = (Shear_Force_1*1000)/((Length_of_Shorter_Side_of_Footing*1000)*Provided_Effective_Depth);
printf("Actual_Shear_Stress = %d N/mm^2 \n", tv1)
tc1 = interp2(tables,tables,tables,Fck,pt);
printf("Permissible_Shear_Stress = %d N.mm^2 \n", tc1)

if (tv1<tc1)
  disp("Hence Safe")
elseif(tv1>tc1)
disp("Hence Unsafe and Shear Reinforcement will be necessary")
Asv1 = no_of_legged*((pi/4)*(dia_stirrups*dia_stirrups));
printf("Area_of_stirrups_reinforcement = %d mm^2 \n", Asv1)
Spacing_of_stirrups = round((0.87*Fy*Asv1*Provided_Effective_Depth)/(Shear_Force_1*1000)/10)*10;
printf("Spacing_of_stirrups= %d mm c/c \n", Spacing_of_stirrups)
endif

disp("\n")
disp("Check for One Way Shear Diagonal Tension For Column A and Column B")
Shear_Force_2 = 1.5* Shear_Force_at_First_Point_of_Contraflexure;
printf("Shear_Force = %d KN \n", Shear_Force_2)

tv2 = (Shear_Force_2*1000)/(Length_of_Shorter_Side_of_Footing*1000*Provided_Effective_Depth);
printf("Actual_Shear_Stress = %d N/mm^2 \n", tv2)

tc2 = interp2(tables,tables,tables,Fck,pt_1);
printf("Permissible_Shear_Stress = %d N/mm^2 \n", tc2)

if (tv2<tc2)
  disp("Hence Safe")
elseif(tv2>tc2)
disp("Hence Unsafe and Shear Reinforcement will be necessary")
Asv2 = no_of_legged*((pi/4)*(dia_of_stirrups*dia_of_stirrups));
printf("Area_of_stirrups_reinforcement = %d mm^2 \n", Asv2)
Spacing_of_stirrups = round((0.87*Fy*Asv2*Provided_Effective_Depth)/(Shear_Force_2*1000)/10)*10;
printf("Spacing_of_stirrups= %d mm c/c \n", Spacing_of_stirrups)
endif

disp("\n")
disp("Check for One Way Shear Diagonal Tension For Column A")
d2 = Width1_of_Column_A/2;
Shear_Force_3 = 1.5*(Shear_Force_just_right_of_A-(Pressure_per_meter_length*d2));
printf("Spacing_of_stirrups= %d KN \n", Shear_Force_3)

tv3 = (Shear_Force_3*1000)/(Length_of_Shorter_Side_of_Footing*1000*Provided_Effective_Depth);
printf("Actual_Shear_Stress = %d N/mm^2 \n", tv3)

tc3 = interp2(tables,tables,tables,Fck,pt_1);
printf("Permissible_Shear_Stress = %d N/mm^2 \n", tc3)

if (tv3<tc3)
  disp("Hence Safe")
elseif(tv3>tc3)
disp("Hence Unsafe and Shear Reinforcement will be necessary")
Asv3 = no_of_legged*((pi/4)*(dia_stirrups*dia_stirrups));
printf("Area_of_stirrups_reinforcement = %d mm^2 \n", Asv3)
Spacing_of_stirrups = round((2.175*Asv3*Fy)/(Length_of_Shorter_Side_of_Footing*1000)/10)*10;
printf("Spacing_of_stirrups= %d mm c/c \n", Spacing_of_stirrups)
endif


disp("\n")
% Transverse Reinforcement 
disp("Transverse Reinforcement for Column A")
Projection_a_beyond_the_face_of_column_A = 0.5*(Length_of_Shorter_Side_of_Footing*1000-Width1_of_Column_A*1000);
printf("Projection_a_beyond_the_face_of_column_A = %d mm \n", Projection_a_beyond_the_face_of_column_A)
Width_B1_of_bending_strip = (Width1_of_Column_A*1000 + 2*Provided_Effective_Depth);
printf("Width_B1_of_bending_strip = %d mm \n", Width_B1_of_bending_strip)
Net_upward_pressure_p0 = (Load_on_Column_A)/(Length_of_Shorter_Side_of_Footing*(Width_B1_of_bending_strip/1000));
printf("Net_upward_pressure = %d KN/m^2 \n", Net_upward_pressure_p0)

Moment1 = Net_upward_pressure_p0*(((Projection_a_beyond_the_face_of_column_A/1000)*(Projection_a_beyond_the_face_of_column_A/1000))/2);
printf("Moment = %d KNm \n", Moment1)
Factored_Moment1 = 1.5*Moment1;
printf("Factored_Moment = %d KNm \n", Factored_Moment1)

depth_1 = sqrt((Factored_Moment1*1000000)/(Ru*1000));
Eff_dep_1  = Required_Effective_Depth - dia_12;
Transverse_Ast1 = (((0.5*Fck)/Fy)*(1-sqrt(1-(4.6*Factored_Moment1*1000000)/(Fck*1000*(Eff_dep_1*Eff_dep_1)))))*1000*Eff_dep_1;
printf("Area_of_Steel = %d mm^2 \n", Transverse_Ast1)
Area_of_1_bar = ceil((pi/4)*(dia_12*dia_12));
printf("Area_of_one_bar = %d mm^2 \n", Area_of_1_bar)
Spacing_1 =  ceil((1000*Area_of_1_bar)/(Transverse_Ast1));
printf("Spacing_btw_bars = %d mm c/c \n", Spacing_1)
tbd = interp1 (table(:,1),table(:,2),Fck_for_column);
Development_Lenght_1 = (Fy*dia_12)/(4*tbd*1.6);

printf("Development_Lenght = %d mm \n", Development_Lenght_1 )
Lenght_of_bar_available_1 = (Projection_a_beyond_the_face_of_column_A-Clear_Cover);
printf("Length_of_bar_available = %d mm \n", Lenght_of_bar_available_1)

if(Development_Lenght_1<Lenght_of_bar_available_1)
disp("Hence Safe in Development Lenght")
elseif
disp("Hence not Safe in Development Lenght")
endif
  

disp("\n")
disp("Transverse Reinforcement for Column B")
Projection_a_beyond_the_face_of_column_B = 0.5*(Length_of_Shorter_Side_of_Footing*1000-Width1_of_Column_B*1000);
printf("Projection_a_beyond_the_face_of_column_B = %d mm \n", Projection_a_beyond_the_face_of_column_B)
Width_B2_of_bending_strip = (Width1_of_Column_B*1000 + 2*Provided_Effective_Depth);
printf("Width_B2_of_bending_strip = %d mm \n", Width_B2_of_bending_strip)

Net_upward_pressure_p1 = (Load_on_Column_B)/(Length_of_Shorter_Side_of_Footing*(Width_B2_of_bending_strip/1000));
printf("Net_upward_pressure = %d KN/m^2 \n", Net_upward_pressure_p1)
Moment2 = Net_upward_pressure_p1*(((Projection_a_beyond_the_face_of_column_B/1000)*(Projection_a_beyond_the_face_of_column_B/1000))/2);
printf("Moment = %d KNm \n", Moment2)
Factored_Moment2 = 1.5*Moment2;
printf("Factored_Moment = %d KNm \n", Factored_Moment2)

Transverse_Ast2 = (((0.5*Fck)/Fy)*(1-sqrt(1-(4.6*Factored_Moment2*1000000)/(Fck*1000*(Eff_dep_1*Eff_dep_1)))))*1000*Eff_dep_1;
printf("Area_of_Steel = %d mm^2 \n", Transverse_Ast2)
Area_of_1_bar = ceil((pi/4)*(dia_12*dia_12));
printf("Area_of_one_bar = %d mm^2 \n", Area_of_1_bar)
Spacing_2 =  ceil((1000*Area_of_1_bar)/(Transverse_Ast2));
printf("Spacing_btw_bars = %d mm c/c \n", Spacing_2)
tbd = interp1 (table(:,1),table(:,2),Fck_for_column);
Development_Lenght_2 = (Fy*dia_12)/(4*tbd*1.6);

printf("Development_Lenght = %d mm \n", Development_Lenght_2)
Lenght_of_bar_available_2 = (Projection_a_beyond_the_face_of_column_B-Clear_Cover);
printf("Length_of_bar_available = %d mm \n", Lenght_of_bar_available_2)

if(Development_Lenght_2<Lenght_of_bar_available_2)
disp("Hence Safe in Development Lenght")
elseif(Development_Lenght_2>Lenght_of_bar_available_2)
disp("Hence not Safe in Development Lenght")
endif

disp("\n")
disp("Transverse Reinforcement for Rest of the Footing")
Area_of_steel_for_transverse_reinforcement_for_rest_of_footing = (0.12/100)*(1000*Overall_Depth);
printf("Area_of_steel_for_transverse_reinforcement_for_rest_of_footing = %d mm^2 \n", Area_of_steel_for_transverse_reinforcement_for_rest_of_footing)

Spacing_of_reinforcement_for_rest_of_the_footing = ceil((1000*Area_of_1_bar)/(Area_of_steel_for_transverse_reinforcement_for_rest_of_footing)/10)*10;
printf("Spacing_of_reinforcement_bars_for_rest_of_the_footing = %d mm c/c \n", Spacing_of_reinforcement_for_rest_of_the_footing)

