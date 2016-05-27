%% Make screenshot of the current figure
function make_snapshot(scrsht_name)
    snapshot = getframe(gcf);
    imgg = frame2im(snapshot);
    imwrite(imgg,scrsht_name);
end