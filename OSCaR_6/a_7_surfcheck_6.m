function a_7_surfcheck_6(const, drow, dpile, surface)

surf(surface.xq,surface.yq,surface.zog);
shading interp
view(2)
hold on
scatter3(drow.npx,drow.npy,surface.F_og(drow.npx,drow.npy),'red')
scatter3(drow.npx,drow.spy,surface.F_og(drow.npx,drow.spy),'black')
daspect([1 1 1])
foo=1;