function plotMatrix(id, A)
  % Visualización matrices
  hf = figure(id)
  A(A == 0) = NaN;
  imagesc(A);
  colorbar;
  axis equal tight;
  colormap(jet);
  title(['Matriz de rigidez global (', num2str(id), ')' ] );
  print (hf, ['/home/santi/UNLP/SistDinamicos/Mod2/parte2/presentacion/ensamble/plotMatrix_', num2str(id), '.jpg'], "-djpeg");
end