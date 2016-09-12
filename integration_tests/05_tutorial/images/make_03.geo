General.Trackball = 0;
General.RotationX = 270;
General.RotationZ = 270;
General.TranslationX = 0;
General.TranslationY = 30;

General.GraphicsWidth = General.MenuWidth + 1400;
General.GraphicsHeight = 1000;
General.ScaleX = 3;
General.ScaleY = 3;
General.ScaleZ = 3;

General.SmallAxes = 0;

Merge "../ref_out/03_tunnel/flow.msh";
View[0].Visible = 1;
View[0].RangeType = 2;
View[0].CustomMin = 0;
View[0].CustomMax = View[0].Max;
View[1].Visible = 0;
View[2].Visible = 1;
View[2].NormalRaise = 1;
View[3].Visible = 1;
View[3].ScaleType = 2;
View[3].ArrowSizeMax = 40;
View[3].NormalRaise = 1e5;
Mesh.SurfaceEdges = 1;
Mesh.SurfaceFaces = 0;

Draw;
Print Sprintf("03_flow.pdf");



Merge "../ref_out/03_tunnel/transport.msh";
View[0].Visible = 0;
View[2].Visible = 0;
View[3].Visible = 0;
View[4].Visible = 1;
View[4].TimeStep = View[4].NbTimeStep-1;
View[4].RangeType = 2;
View[4].CustomMin = -10.7;
View[4].CustomMax = -10.3;
View[4].SaturateValues = 1;
Mesh.SurfaceEdges = 1;
Mesh.SurfaceFaces = 0;

Draw;
Print Sprintf("03_transport.pdf");
Exit;