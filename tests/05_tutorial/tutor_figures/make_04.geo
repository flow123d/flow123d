General.Trackball = 0;
General.RotationX = -90;
General.RotationZ = 0;
General.TranslationX = -35;
General.TranslationY = 0;

General.GraphicsWidth = General.MenuWidth + 1000;
General.GraphicsHeight = 1000;
General.ScaleX = 3.5;
General.ScaleY = 3.5;
General.ScaleZ = 3.5;

General.SmallAxes = 0;
General.Axes = 1;

//Merge "../04_mesh.msh";
Mesh.Lines = 0;
Mesh.SurfaceEdges = 1;
Mesh.SurfaceFaces = 0;
Mesh.ColorCarousel = 0; // color by physical labels

Merge "../ref_out/04_frac_diffusion/flow.msh";
View[0].DrawLines = 1;
View[0].DrawTriangles = 0;
View[0].LineWidth = 5;
View[0].ColormapNumber = 0;
View[0].ShowScale = 0;
View[1].Visible = 0;
View[2].Visible = 0;
View[3].Visible = 0;


Draw;
Print Sprintf("04_mesh.pdf");


Exit;

