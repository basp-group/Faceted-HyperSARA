function get_residual(src_filename, dest_filename)

    m = matfile(src_filename);
    fitswrite(m.res, strcat(dest_filename, ".fits"));

end