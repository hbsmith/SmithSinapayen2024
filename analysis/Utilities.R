# This file contains all functions we ever implemented for data analysis.

library(rjson)
library(kneedle)
library(ecodist)

# general utilities
alive_value = 1 

# shape = "1col" or "2col"
prepare_figure_pdf = function(figure_folder, file_name, param_names,  shape){
  figure_path = paste0(figure_folder, file_name, "_", param_names, ".pdf")
  print(paste("Saving figure at", figure_path))
  
  if (shape=="1col"){
    w = 5
    h = 4
  } else {
    w = 8
    h = 3  
  }
  
  pdf(file = figure_path,   
      width = w,
      height = h)
}


# read model, planer data and parameters
# remove NAs
read_sim_data = function(base_path, sub_folder){
  planets_path = paste(base_path, sub_folder, "df_planets.csv.gz", sep="")
  life_path = paste(base_path, sub_folder, "df_lifes.csv.gz", sep="")
  parameters_path = paste(base_path, sub_folder, "parameters.json", sep="")
  
  parameters = fromJSON(file=parameters_path)
  
  zz = gzfile(planets_path,'rt')  
  df_planets_sweep = read.csv(zz, stringsAsFactors=F, sep=";")
  close(zz)
  
  zz = gzfile(life_path,'rt')  
  df_lifes_sweep = read.csv(zz, stringsAsFactors=F, sep=";")
  close(zz)
  beep()
  
  model_path = paste(base_path, sub_folder, "df_model.csv.gz", sep="")
  zz = gzfile(model_path,'rt')  
  df_model_sweep = read.csv(zz, stringsAsFactors=F, sep=";")
  close(zz)
  # replace NAN with 0
  df_model_sweep[is.na(df_model_sweep)] = 0
  
  l = list(parameters, df_planets_sweep, df_lifes_sweep, df_model_sweep)
  return(l)
}


# calculate physical distance
# returns 2 matrices
# m0 is not normalized
# m1 is normalized
physical_distance = function(planets, space_dim=3){
  if (space_dim == 2){
    coord_names = c("x","y")
  } else {
    coord_names = c("x","y","z")
  }
  
  coordinates = data.matrix(planets[coord_names])
  distances = full(distance(coordinates))
  distances_norm = distances/max(distances)
  return(list(distances, distances_norm))
}

# returns matrices with all composition distances
# m0 is not normalized
# m1 is normalized
composition_distance = function(planets, comp_dim=10){
  # calculate distance between all composition vectors
  comp_names = paste("comp_",c(1:comp_dim), sep="")
  compositions = data.matrix(planets[comp_names])
  comp_dist = full(distance(compositions))
  distances_norm = comp_dist/max(comp_dist)
  return(list(comp_dist, distances_norm))
}


plot_comp_dist = function(planets, title){
  # color plot
  pal <- colorRampPalette(rev(brewer.pal(n = 9, name = "Blues")))
  col_order = findInterval(planets$comp_dist, sort(planets$comp_dist))
  plot(planets$x, planets$y, main=title, xlab="x", ylab="y", pch=21, col="black", 
       # bg="red",
       bg=pal(nrow(planets))[col_order])
  points(start_planet$x, start_planet$y, pch=19)
}

# tests to run to check data integrity

# Plot composition distributions
plot_distributions = function(planets_t0, planets_t){
  n_terraformed = length(subset(planets_t, alive==alive_value)[,1])
  title = "Distribution change should reflect number of terraformed planets"
  hist(planets_t0$comp_distance,
       xlab = "Composition Distance",
       main = title)
  title = paste(n_terraformed, " total terraformed planets", sep="")
  hist(planets_t$comp_distance,
       xlab = "Composition Distance",
       main = title)
}


plot_trajectories = function(planets_t0, planets_t){
  # use when saving
  # output_path = "./figures/"
  # png(file.path(output_path, "before.png"), units="mm", width=200, height=150, res=300)
  plot_comp_dist(planets_t0, "the color change should make sense")
  # dev.off()
  plot_comp_dist(planets_t, "Blue = high composition distance")
  points(df_lifes$x, df_lifes$y, type = "l", col="gray")
}


plot_distances = function(planets_t0, planets_t, title="Normalized distances"){
  plot(planets_t0$distance, planets_t0$comp_distance, col="gray",# type = "l", 
       xlab = "Physical distance, normalised", ylab = "Composition distance",
       main = title
  )
  plot(x=planets_t$distance, y=planets_t$comp_distance, col="blue",# type = "l", 
       xlab = "Physical distance, normalised", ylab = "Composition distance",
       main = title
  )
}

# counts number terraformed planets at timestep t
# t i the timestep
# type_var: the value of alive (eg "True" or True or 1)
count_live = function(t, df){
  sub = subset(df, step==t)
  alive_count = sum(sub$alive==alive_value)
  return(alive_count)
}

# p_distance_mat is 
get_mantel_fast = function(steps, planets, p_distance_mat, verbose=TRUE){
  mantels = data.frame(matrix(nrow = 0, ncol = 5))
  max_step = tail(planets$step, 1)
  
  b_planets_t = subset(planets, step == steps[1])
  n_planets = max(b_planets_t$id)
  terraformed = subset(b_planets_t, alive == alive_value)
  all_terraformed = terraformed
  
  c_distance_mat = lower(composition_distance(b_planets_t)[[1]])
  
  for (t in steps){
    if(verbose){
      print(paste0(t, " out of ", max_step))
    }
    
    if (t != steps[1]){
      # subset for t
      b_planets_t = subset(planets, step == t)
      # find which planet changed
      all_terraformed = subset(b_planets_t, alive == alive_value)
      new_terraformed = subset(all_terraformed, !id %in% terraformed$id)
      
      # calculate comp distance
      new_terraformed_c = new_terraformed[comp_names]
      comps = b_planets_t[comp_names]
      
      for (i in c(1:length(new_terraformed$id))){
        nc = t(new_terraformed_c[i,])
        temp = sweep(comps, 2, nc, `-`)
        temp = temp**2
        temp = apply(temp, 1, sum)
        c_dist = sqrt(temp)
        
        # update the composition distance matrix
        temp = full(c_distance_mat)
        temp[new_terraformed$id[i], ] = c_dist
        temp[ ,new_terraformed$id[i]] = c_dist
        c_distance_mat = lower(temp)
      }
      
    }
    
    # Mantel test: is the geographic distance between planets
    # related to the (difference in) composition between planets
    # error when all compositions are the same (mantel coeff = 1)
    tryCatch(
      expr = {
        m = mantel(c_distance_mat ~ p_distance_mat, nperm=n_perm)
      },
      error = function(e){ 
        print(e)
      },
      warning = function(w){
        print(w)
      },
      finally = {
      }
    )
    
    # m = mantel(c_distance_mat ~ p_distance_mat, nperm=n_perm)
    mantels = rbind(mantels, c(t,m[1:4]))
    
    terraformed = all_terraformed
  }
  
  colnames(mantels) = c("t", "mantelr", "pval1", "pval2", "pval3")
  return(mantels)
}


terraformed_ratio = function(df_planets){
  # extract relevant time steps
  max_steps = max(df_planets$step)
  unique_steps = unique(df_planets$step)
  alive_count = sapply(unique_steps, count_live, df=df_planets)
  alive_count = alive_count/length(unique(df_planets$id))
  alive_count = data.frame(unique_steps, alive_count)
  colnames(alive_count) = c("t", "alive_ratio")
  return(alive_count)
}

# Fig 1 utilities

# Plot different subdivisions of model for visual clustering
# divisions: c(name of parameters to plot)
plot_divisions = function(df_model_sweep, parameters, div_colnames){
  
  plot(df_model_sweep$n_living_planets/1000, df_model_sweep$mantel_corr_coeff, main="All data")
  
  max_y = 0.4
  
  for(c_name in div_colnames){
    
    plot(1, type = "n", main = c_name,
         xlab = "% terraformed", ylab = "mantel",
         xlim=c(0,1), ylim=c(0, max_y) 
    )
    
    col = 1
    for (n_div in parameters[[c_name]]){
      df_model = subset(df_model_sweep, df_model_sweep[c_name]==n_div)
      points(df_model$n_living_planets/1000, df_model$mantel_corr_coeff, col=col)
      col = col+1
    }
  }
}

# subdivision is one of c("nool", "n_idxs_to_keep_from_destination", "allowed_diff", "r", "mutation_rate")
plot_subdivision = function(subdivision, df_model_sweep, parameters, divisions){
  
  # the chosen subdivision
  for (div in parameters[[subdivision]]){
    print(paste(subdivision, div))
    
    df_model_div = subset(df_model_sweep,  df_model_sweep[subdivision]==div)
    divisions_reduced = divisions[ !divisions == subdivision]
    
    plot_divisions(df_model_div, parameters, divisions_reduced)
  }
}


get_tile_centers = function(param){
  l = length(param)
  # this tile is
  # middle of the space between the previous tile and the next
  # + previous tile's end as offset
  # start of first tile
  start = 0 
  bottom_half = (param[1] - start)
  top_half = (param[2] - param[1])/2
  a = (bottom_half + top_half)/2
  
  tile_center = c(a)
  tile_size = c(bottom_half + top_half)
  for(i in 2:l-1){
    bottom_half = (param[i] - param[i-1])/2
    top_half = (param[i+1] - param[i])/2
    t_size = bottom_half + top_half
    # add offset and find middle
    a = param[i] - bottom_half + t_size/2 
    tile_center = c(tile_center, a)
    tile_size = c(tile_size, t_size)
  }
  
  i = l
  bottom_half = (param[i] - param[i-1])/2
  # no top half
  top_half = 0
  t_size = bottom_half + top_half
  # add offset and find middle
  a = param[i] - bottom_half + t_size/2
  tile_center = c(tile_center, a)
  tile_size = c(tile_size, t_size)
  
  df = data.frame(tile_center, tile_size)
  
  return(df)
}


# selction_function = one of "r" or "diff"
# n_ool = 1 or 10
make_heatmap_data = function(my_df_model, my_selection_function, my_parameters, my_n_ool, use_pvalue = FALSE){
  heatmap_data = data.frame(matrix(nrow = 0, ncol = 3)) 
  df_model_ool = subset(my_df_model,  nool==my_n_ool)
  compsize = my_parameters$compsize
  
  if (my_selection_function == "diff"){
    columns = c("n_id" ,"n_diff", "max_mr", "min_mr") 
    
    for (n_id in my_parameters$n_idxs_to_keep_from_destination){
      df_model_nid = subset(df_model_ool, n_idxs_to_keep_from_destination==n_id)
      
      for (n_diff in my_parameters$allowed_diff){
        df_model = subset(df_model_nid, allowed_diff==n_diff)
        if(use_pvalue){
          max_mr = max(df_model$mantel_p_value)
          min_mr = min(df_model$mantel_p_value)
        } else {
          max_mr = max(df_model$mantel_corr_coeff)
          min_mr = min(df_model$mantel_corr_coeff)
        }
        data = c(n_id, n_diff, max_mr, min_mr)
        heatmap_data = rbind(heatmap_data, data)
      }
    }
  } else {
    columns = c("n_id","n_r","max_mr", "min_mr")  
    
    for (n_id in my_parameters$n_idxs_to_keep_from_destination){
      df_model_nid = subset(df_model_ool, n_idxs_to_keep_from_destination==n_id)
      
      for (n_r in my_parameters$r){
        df_model = subset(df_model_nid, r==n_r)
        max_mr = max(df_model$mantel_corr_coeff)
        min_mr = min(df_model$mantel_corr_coeff)
        data = c(n_id, n_r, max_mr, min_mr)
        heatmap_data = rbind(heatmap_data, data)
      }
    }
  }
  
  colnames(heatmap_data) = columns
  return(heatmap_data)
}


# calculates cell sizes
proportional_heatmap_data = function(my_heatmap_data, my_parameters, my_selection_function){
  
  a = get_tile_centers(my_parameters$n_idxs_to_keep_from_destination)
  
  if (my_selection_function == "diff"){
    a = data.frame(rep(a$tile_center, each=length(my_parameters$allowed_diff)), rep(a$tile_size, each=length(my_parameters$allowed_diff)))
    colnames(a) = c("tile_center", "tile_size")
    
    b = get_tile_centers(my_parameters$allowed_diff)
    l = length(my_parameters$n_idxs_to_keep_from_destination)
    b = data.frame(rep(b$tile_center, l), rep(b$tile_size, l))
    colnames(b) = c("tile_center", "tile_size")
    
  } else {
    a = data.frame(rep(a$tile_center, each=length(my_parameters$r)), rep(a$tile_size, each=length(my_parameters$r)))
    colnames(a) = c("tile_center", "tile_size")
    
    b = get_tile_centers(my_parameters$r)
    l = length(my_parameters$n_idxs_to_keep_from_destination)
    b = data.frame(rep(b$tile_center, l), rep(b$tile_size, l))
    colnames(b) = c("tile_center", "tile_size")
  }
  
  my_heatmap_data$tile_x = a$tile_center
  my_heatmap_data$tile_y = b$tile_center
  my_heatmap_data$tile_w = a$tile_size
  my_heatmap_data$tile_h = b$tile_size
  return(my_heatmap_data)
}


# heatmap_data parameter must be already augmented with tile size
draw_proportional_heatmap = function(my_heatmap_data, my_parameters, my_selection_function, legend=TRUE,
                                     start_low=TRUE, use_pvalue = FALSE){
  if (use_pvalue){
    legend_name = "Max p-value"
  } else {
    legend_name = "Max Mantel"
  }
  
  # for color scale
  if (use_pvalue){
    rng = c(0, max(my_heatmap_data$max_mr))
  } else {
    rng = c(0, 0.5)
  }
  comp_percent = round(100*my_parameters$n_idxs_to_keep_from_destination/my_parameters$compsize, digits = 0)
  labels_x = paste0(comp_percent, "%")
  
  if (my_selection_function == "diff"){
    max_diff = round(sqrt(my_parameters$compsize), digits=1)
    diff_percent = round(100*my_parameters$allowed_diff/max_diff, digits = 0)
    labels_y = paste0(diff_percent, "%")

    ggp <- ggplot(my_heatmap_data, show.legend = legend,
                  aes(x=tile_x, y=tile_y,
                      width = tile_w,
                      height = tile_h
                  )
    )
    if(use_pvalue){
      ggp = ggp + geom_tile(aes(fill = min_mr), alpha = 1)
    } else {
      if (start_low){
        ggp = ggp + geom_tile(aes(fill = max_mr), alpha = 1) 
      } else {
        ggp = ggp + geom_tile(aes(fill = min_mr), alpha = 1)
      }
    }
    
    ggp = ggp +                    
      scale_x_continuous(breaks = my_parameters$n_idxs_to_keep_from_destination,
                         minor_breaks = NULL,
                         labels = labels_x) + # waiver() 
      scale_y_continuous(breaks = my_parameters$allowed_diff,
                         minor_breaks = NULL,
                         labels = labels_y
      ) +
      scale_fill_gradient2(low = "white", mid = "#619CFF", high="black", midpoint=mean(rng),    #same midpoint for plots (mean of the range)
                           breaks=seq(rng[1], rng[2], 0.25), #breaks in the scale bar
                           limits=rng, name=legend_name
      )
    if(legend){
      ggp = ggp + xlab(label = "Retained planetary composition") + 
        ylab(label = "Maximum allowed composition difference")
    } else {
      ggp = ggp + xlab(label = element_blank()) +  ylab(label = element_blank())
    }
  } else {
    # what is the real max value for r?
    max_r = max(my_parameters$r)
    r_percent = round(100*my_parameters$r/max_r, digits = 0)
    labels_y = paste0(r_percent, "%")
    
    ggp <- ggplot(my_heatmap_data, 
                  aes(x=tile_x, y=tile_y,
                      width = tile_w,
                      height = tile_h
                  )
    )
    if (start_low){
      ggp = ggp + geom_tile(aes(fill = max_mr), alpha = 1) 
    } else {
      ggp = ggp + geom_tile(aes(fill = min_mr), alpha = 1)
    }
    ggp = ggp +  
      scale_x_continuous(breaks = my_parameters$n_idxs_to_keep_from_destination,
                         minor_breaks = NULL,
                         labels = labels_x ) +
      scale_y_continuous(breaks = my_parameters$r,
                         minor_breaks = NULL,
                         labels = labels_y
      ) +
      scale_fill_gradient2(low = "white", mid = "#619CFF", high="black", midpoint=mean(rng),    #same midpoint for plots (mean of the range)
                           breaks=seq(rng[1], rng[2], 0.25), #breaks in the scale bar
                           limits=rng)
    if(legend){
      ggp = ggp + xlab(label = "Retained planetary composition") +
        ylab(label = "allowed radius")
    } else {
      ggp = ggp + xlab(label = element_blank()) +  ylab(label = element_blank())
    }
  }
  
  return(ggp)
}


# returns the heatmap with highlighted cell
# selection_function: one of c( "r", "diff")
heatmap_with_highlight = function(my_df_model, my_selection_function, my_chosen_dfs, my_parameters, my_n_ool, my_colors, legend=TRUE){
  
  heatmap_data = make_heatmap_data(my_df_model, my_selection_function, my_parameters, my_n_ool)
  heatmap_data = proportional_heatmap_data(heatmap_data, my_parameters, my_selection_function)
  ggp = draw_proportional_heatmap(heatmap_data, my_parameters, my_selection_function, legend = legend)
  #ggp
  
  # add rectangles
  #  sizes
  n_col_name = paste0("n_",my_selection_function)
  tile_sizes = data.frame(x = unique(heatmap_data$tile_x), y=unique(heatmap_data$tile_y),
                          w = unique(heatmap_data$tile_w), h=unique(heatmap_data$tile_h),
                          id = unique(heatmap_data$n_id), selection_function = unique(heatmap_data[[n_col_name]]))
  
  sub_df = subset(my_chosen_dfs, ool == my_n_ool & selection==my_selection_function)
  l = length(sub_df$selection)
  sub_df$x = rep(0, l)
  sub_df$y = rep(0, l)
  sub_df$h = rep(0, l)
  sub_df$w = rep(0, l)
  
  cell_labels = rep("", length(heatmap_data$max_mr))
  if (my_selection_function == "diff"){
    col_vals = my_parameters$allowed_diff
  } else {
    col_vals = my_parameters$r
  }
  row_vals = my_parameters$n_idxs_to_keep_from_destination
  
  for (i in c(1:l)){
    sub_df$w[i] = 0+subset(tile_sizes, id==sub_df$id[i])$w
    sub_df$h[i] = 0+subset(tile_sizes, selection_function==sub_df$select_value[i])$h
    # we want starting points, not centers
    sub_df$x[i] = subset(tile_sizes, id==sub_df$id[i])$x - sub_df$w[i]/2
    sub_df$y[i] = subset(tile_sizes, selection_function==sub_df$select_value[i])$y - sub_df$h[i]/2
    # label only for these cells
    if (my_selection_function == "diff"){
      max_mr = subset(heatmap_data, n_id == sub_df$id[i] & n_diff == sub_df$select_value[i])$max_mr
    } else {
      max_mr = subset(heatmap_data, n_id == sub_df$id[i] & n_r == sub_df$select_value[i])$max_mr
    }
    lab_index = (which(row_vals==sub_df$id[i]) - 1) * length(col_vals)  + which(col_vals==sub_df$select_value[i])
    cell_labels[lab_index] = round(max_mr, 2) #paste0("",round(max_mr, 2))
  }
  
  ggp_high = ggp +
    annotate("rect", xmin=sub_df$x, xmax=sub_df$x+sub_df$w,
             ymin=sub_df$y, ymax=sub_df$y + sub_df$h,
             fill=NA, color=my_colors[(sub_df$curve_type+1)], size = 2) +
    # geom_text(aes(label = cell_labels, colour="red"), size=3.5) # add labels
    geom_text(aes(label = cell_labels, colour = ifelse(cell_labels<0.25, "black", "white"))) +
    scale_colour_manual(values=c("white"="white", "black"="black"))+
    theme(legend.position="none")
  
  # ggp_high
  
  return(ggp_high)
}

plot_selected_mantels = function(my_df_model, my_chosen_dfs, my_selection_function){
  
  n_cols = length(colnames(my_df_model)) + 1
  sub_df = data.frame(matrix(0, nrow = 0, ncol = n_cols))
  
  # iterate on columns
  for (i in c(1:length(my_chosen_dfs$ool))){
    if(my_selection_function == "diff") {
      if (my_chosen_dfs$selection[i] == "diff"){
        temp = subset(my_df_model, nool == my_chosen_dfs$ool[i]
                      & allowed_diff == my_chosen_dfs$select_value[i]
                      & mutation_rate == my_chosen_dfs$mutation[i]
                      & n_idxs_to_keep_from_destination == my_chosen_dfs$id[i]
        )
        
        # you might choose the color here
        temp$curve_type = my_chosen_dfs$curve_type[i]
        #
        sub_df = rbind(sub_df, temp)
      }
    } else {
      if (my_chosen_dfs$selection[i] == "r"){
        temp = subset(my_df_model, nool == my_chosen_dfs$ool[i]
                      & r == my_chosen_dfs$select_value[i]
                      & mutation_rate == my_chosen_dfs$mutation[i]
                      & n_idxs_to_keep_from_destination == my_chosen_dfs$id[i]
        )
        
        # you might choose the color here
        temp$curve_type = my_chosen_dfs$curve_type[i]
        #
        sub_df = rbind(sub_df, temp)
      }
    }
  }
  
  colnames(sub_df) = c(colnames(my_df_model), "curve_type")
  
  # how about plot everything in grey then this on top
  plot_all = ggplot(my_df_model, aes(n_living_planets/1000, mantel_corr_coeff)) +
    geom_point(color = "grey") + labs(x = "Alive ratio", y = "Mantel ratio")
  
  plot_selected = plot_all + geom_point(data = sub_df, color = colors[(sub_df$curve_type+1)])
  return(plot_selected)
}


# New compact design
make_heatmap_figure = function(my_selection_function, my_heatmap_1, my_heatmap_10, my_mantels_plot, my_save=FALSE, my_filename = ""){
  # https://github.com/cxli233/Online_R_learning/blob/master/Quick_data_vis/Lessons/08_Intro_plot_assembly.md
  
  design <- c(
    "A
    B"
  )
  
  # TOOD: save somewhere else?
  if (my_save) {
    pdf(file = my_filename,   
        width = 15,
        height = 9)
    # png(paste0(my_base_path, my_filename), units="px", width=4000, height=2000, res=300)
  }
  
  # combined heatmaps
  sub_design = c(
    "AB"
  )
  if(my_selection_function == "diff"){
    y_lab = "Maximum allowed composition difference"
  } else{
    y_lab = "Maximum allowed distance (normalized)"
  }
  Panel_1 = wrap_plots(my_heatmap_10 + theme(legend.position = "none")  + 
                         labs(tag = "10 OoL")  + 
                         xlab("Retained planetary composition") + ylab(y_lab),
                       my_heatmap_1 + labs(tag = "1 OoL"),
                       design = sub_design
  )
  
  fig = wrap_plots(Panel_1,
                   my_mantels_plot,
                   design = design
  )
  
  if (my_save){
    print(fig)
    dev.off()
    print(paste("Saved at", my_filename))
  }
  
  return(fig)
}

# Calculate test accuracy
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4614595/#:~:text=Accuracy%3A%20The%20accuracy%20of%20a,TP%20%2B%20TN%20%2B%20FP%20%2B%20FN
# stats is the result of xxx function or
# names(stats) = c("true_positives", "false_positives", "true_negatives", "false_negatives")
test_eval = function(stats){
  tp = stats["true_positives"]
  tn = stats["true_negatives"]
  fp = stats["false_positives"]
  fn = stats["false_negatives"]
  
  accuracy = (tp + tn) / sum(stats)
  sensitivity = tp / (tp + fn)
  specificity = tn / (tn + fp)
  
  df = data.frame(accuracy, sensitivity, specificity)
  rownames(df) <- NULL
  return(df)
}


get_planets_betweeen = function(df_planets, min_x, max_x, coord_name){
  cell = subset(df_planets, df_planets[coord_name]>=min_x & df_planets[coord_name]<=max_x)
  return(cell)
}


## Plot all curves to pick representative ones
plot_runs_individually = function(parameters, df_model_sweep, selection_func){
  
  ids = parameters$n_idxs_to_keep_from_destination
  mutations = parameters$mutation_rate
  diffs = parameters$allowed_diff
  ools = parameters$nool
  rs = parameters$r
  
  if (selection_func == "r"){
    param_order = "ool r mut id"
    sf = rs
  } else {
    param_order = "ool diff mut id"
    sf = diffs
  }
  
  for (n_ool in ools){
    df_model_ool = subset(df_model_sweep,  nool==n_ool)
    
    for (n_diff in sf){
      if (selection_func == "r"){
        df_model_diff = subset(df_model_ool,  r==n_diff)
      } else {
        df_model_diff = subset(df_model_ool,  allowed_diff==n_diff)
      }
      for (n_mutation in mutations){
        df_model_mut = subset(df_model_diff,  mutation_rate==n_mutation)
        for (n_id in ids){
          df_model_id = subset(df_model_mut,  n_idxs_to_keep_from_destination==n_id)
          title = paste(param_order, n_ool, n_diff, n_mutation, n_id)
          
          plot(df_model_id$n_living_planets/1000, df_model_id$mantel_corr_coeff,
               ylim = c(0,0.5), xlim = c(0,1),
               main = title
          )
        }
      }
    }
  }
}


# detected = subset from df_planets or at least df with id column
evaluate_detection = function(detected, my_df_planets){
  alive_planets = subset(my_df_planets, alive == alive_value)
  not_alive_planets = subset(my_df_planets, alive != alive_value)
  not_detected = subset(my_df_planets, !(id %in% detected$id))
  
  # true positives
  temp = subset(detected, id %in% alive_planets$id)
  tp = length(temp$x)
  
  # false positives
  temp = subset(detected, id %in% not_alive_planets$id)
  fp = length(temp$x)
  
  # true negatives
  temp = subset(not_detected, id %in% not_alive_planets$id)
  tn = length(temp$x)
  
  # false negatives
  temp = subset(not_detected, id %in% alive_planets$id)
  fn = length(temp$x)
  
  stats = c(tp, fp, tn, fn)
  names(stats) = c("true_positives", "false_positives", "true_negatives", "false_negatives")
  eval = test_eval(stats)
  results = c(stats, eval)
  
  return(results)
}


select_parameter_set = function(df_planet_sweep, df_life_sweep, n_id, n_diff, n_mutation, selection_func){
  print("selection function")
  print(selection_func)
  print("n_idxs_to_keep_from_destination")
  print(n_id)
  print("allowed_radius")
  n_diff = p_r
  print(n_diff)
  print("mutation_rate")
  n_mutation = p_mutation
  print(n_mutation)
  
  df_planets = subset(df_planet_sweep, n_idxs_to_keep_from_destination==n_id)
  if (selection_func == "r"){
    df_planets = subset(df_planets, r==n_diff)
  } else {
    df_planets = subset(df_planets, allowed_diff==n_diff)
  }
  df_planets = subset(df_planets, mutation_rate==n_mutation)
  df_planets = subset(df_planets, nool==n_ool)
  
  df_lifes = subset(df_life_sweep, n_idxs_to_keep_from_destination==n_id)
  if (selection_func == "r"){
    df_lifes = subset(df_lifes, r==n_diff)
  } else {
    df_lifes = subset(df_lifes, allowed_diff==n_diff)
  }
  df_lifes = subset(df_lifes, mutation_rate==n_mutation)
  df_lifes = subset(df_lifes, nool==n_ool)
  
  df_model = subset(df_model_sweep, n_idxs_to_keep_from_destination==n_id)
  if (selection_func == "r"){
    df_model = subset(df_model, r==n_diff)
  } else {
    df_model = subset(df_model, allowed_diff==n_diff)
  }
  df_model = subset(df_model, mutation_rate==n_mutation)
  df_model = subset(df_model, nool==n_ool)
  
  if (selection_func == "r"){
    param_order = "ool r mut id"
  } else {
    param_order = "ool diff mut id"
  }
  
  title = paste(sub_folder, param_order, n_ool, n_diff, n_mutation, n_id)
  plot(df_model$n_living_planets/1000, df_model$mantel_corr_coeff,
       ylim = c(0,0.5), xlim = c(0,1),
       main = title, type = "l"
  )
  
  return(list(df_planets, df_lifes, df_model))
}


# Find the epsilon parameter for dbscan
# @param comp_distances distance matrix for compositions
epsilon_dbscan = function(comp_distances, db_minPts, total_n_planets, plot_elbow = FALSE, my_sensitivity = 4){
  distances_elbow = c()
  
  for (i in c(1:total_n_planets)){
    # find closest neighbors
    d = comp_distances[i,]
    df = data.frame(d)
    colnames(df) = c("d")
    df$id = c(1:total_n_planets)
    
    df = df[order(df$d), ]
    neighbor = df$d[db_minPts+1]
    
    distances_elbow = c(distances_elbow, neighbor)
  }
  
  d = sort(distances_elbow)
  df = data.frame(d)
  colnames(df) = c("d")
  df$id = c(1:total_n_planets)
  
  knee <- kneedle(df$id, df$d, sensitivity = my_sensitivity, decreasing = FALSE)
  elbow = knee[2]
  
  if (plot_elbow){
    plot(d, type="l", xlab = "Sorted neighbors", ylab = "Composition distance",
         main = paste("Sensitivity = ", my_sensitivity)
           )
    abline(h = elbow, col = "red", lty = 2)
  }
  
  return(elbow)
}


dbscan_clusters = function(db_minPts, comp_distances, total_n_planets, do_plot_elbow = FALSE, sensitivity = 4){
  elbow = epsilon_dbscan(comp_distances, db_minPts, total_n_planets, plot_elbow = do_plot_elbow, my_sensitivity = sensitivity)
  dbscan_cl = dbscan(compositions, eps = elbow, minPts = db_minPts)
  
  return(dbscan_cl)
}
