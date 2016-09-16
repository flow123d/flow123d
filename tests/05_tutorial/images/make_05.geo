General.Trackball = 0;
General.RotationX = -90;
General.RotationZ = 0;
General.TranslationX = -1700;
General.TranslationY = 0;

General.GraphicsWidth = General.MenuWidth + 1000;
General.GraphicsHeight = 1000;
General.ScaleX = 3;
General.ScaleY = 3;
General.ScaleZ = 3;

General.SmallAxes = 0;

Merge "../ref_out/05_heat/flow.msh";
View[0].Visible = 0;
View[1].Visible = 0;
Mesh.SurfaceEdges = 1;
Mesh.SurfaceFaces = 0;

Draw;
Print Sprintf("05_mesh.pdf");


View[0].Visible = 1;
View[1].Visible = 1;
View[1].ScaleType = 2;
Mesh.SurfaceEdges = 0;

Draw;
Print Sprintf("05_flow.pdf");


Merge "../ref_out/05_heat/heat.msh";
View[0].Visible = 0;
View[1].Visible = 0;
View[2].Visible = 1;
View[2].TimeStep = View[2].NbTimeStep-1;
View[2].LineType = 1;
View[2].LineWidth = 10;
View[2].RangeType = 2;
View[2].CustomMin = 283.15;
View[2].CustomMax = 433.15;
View[2].SaturateValues = 1;

Draw;
Print Sprintf("05_transport.pdf");

Exit;