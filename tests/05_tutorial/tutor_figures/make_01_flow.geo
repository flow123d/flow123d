General.Trackball = 0;
General.RotationX = 270;
General.RotationZ = 270;
General.TranslationX = 0;
General.TranslationY = -4;

General.GraphicsWidth = General.MenuWidth + 350;
General.GraphicsHeight = 933;
General.ScaleX = 6.5;
General.ScaleY = 6.5;
General.ScaleZ = 6.5;

General.SmallAxes = 0;

Merge "../ref_out/01_column/flow.msh";
View[0].Visible = 0;
View[1].Visible = 0;
View[2].Visible = 0;
View[3].Visible = 0;
Mesh.SurfaceEdges = 1;
Mesh.SurfaceFaces = 0;
Draw;
Print Sprintf("01_mesh.pdf");


General.TranslationY = -2;
General.GraphicsWidth = General.MenuWidth + 300;
General.GraphicsHeight = 800;
Mesh.SurfaceEdges = 0;
View[2].Visible = 1;
View[3].Visible = 1;
View[3].ArrowSizeMax = 30;
View[3].CustomMax = 1e-08;
View[3].CustomMin = 1e-09;
View[3].RangeType = 2;
View[3].ShowScale = 0;

Draw;
Print Sprintf("01_flow.pdf");
Exit;