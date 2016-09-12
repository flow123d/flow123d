General.Trackball = 0;
General.RotationX = 270;
General.RotationZ = 270;
General.TranslationX = 0;
General.TranslationY = -2;

General.GraphicsWidth = General.MenuWidth + 300;
General.GraphicsHeight = 800;
General.ScaleX = 6.5;
General.ScaleY = 6.5;
General.ScaleZ = 6.5;

General.SmallAxes = 0;

//Mesh.SurfaceEdges = 1;

For num In {1:4}
  Merge "../ref_out/02_column_transport/transport.msh";

  View[0].TimeStep = 5*num;

  Draw;
  Print Sprintf("02_transport_%01g.pdf", num);

  Delete View[0];
EndFor

Exit;