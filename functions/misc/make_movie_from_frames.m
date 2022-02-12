function make_movie_from_frames(folder_path, save_name)
    % save a movie from frames usin ffmpeg
    
    StringArray = {...,
    sprintf("MOVIENAME='%s'",save_name),...
    "echo $MOVIENAME",...
    " ",...
    "ffmpeg -framerate 3 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p -profile:v high -level 4.2 -crf 17 -vf scale=1920:-2 $MOVIENAME.mp4",...
    };
    
    fid = fopen([folder_path,filesep,'MovieGen.sh'], 'wt');
    fprintf(fid, '%s\n', string(StringArray));
    fclose(fid);
    
    disp("Now go to folder_path via terminal, and run 'bash MovieGen.sh'")

end