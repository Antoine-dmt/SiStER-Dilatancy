function [movie_name] = tectoland_movie(movie_name,imin,imax,istep,framerate)

%movie_name = 'REF_cold';




vidObj = VideoWriter(movie_name,'MPEG-4');
vidObj.Quality = 100;
set(vidObj,'FrameRate',framerate)
open(vidObj);


for i = imin:istep:imax
    
    load(num2str(i))
    fastscatter(xm(im>1),ym(im>1),xim(im>1));axis ij; axis equal
    set(gcf, 'Position', get(0, 'Screensize'));
    axis([5e4 10e4 0 2E4])
    colorbar
    caxis([0 1])
    %tectoland_moho(i,Tbdt,1)
    hold off
    pause(.001)
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
    
end


close(vidObj);
