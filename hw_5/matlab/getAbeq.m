function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1)
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    Aeq_start(1,1:n_order+1) = DerivativePolyCoeffient_helper(n_order, 0, 0);
    Aeq_start(2,1:n_order+1) = DerivativePolyCoeffient_helper(n_order, 1, 0);
    Aeq_start(3,1:n_order+1) = DerivativePolyCoeffient_helper(n_order, 2, 0);
    Aeq_start(4,1:n_order+1) = DerivativePolyCoeffient_helper(n_order, 3, 0);

    beq_start = start_cond';

    %#####################################################
    % p,v,a constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end
    %
    %
    %
    %
    Aeq_end(1, n_all_poly-n_order : n_all_poly) = DerivativePolyCoeffient_helper(n_order, 0, ts(end));
    Aeq_end(2, n_all_poly-n_order : n_all_poly) = DerivativePolyCoeffient_helper(n_order, 1, ts(end));
    Aeq_end(3, n_all_poly-n_order : n_all_poly) = DerivativePolyCoeffient_helper(n_order, 2, ts(end));
    Aeq_end(4, n_all_poly-n_order : n_all_poly) = DerivativePolyCoeffient_helper(n_order, 3, ts(end));

    beq_end = end_cond';
    
    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    %
    %
    %
    %
    carry_idx = 1;
    for mid_idx = 2 : n_seg
        Aeq_wp(mid_idx-1, (mid_idx-1)*8+1:(mid_idx)*8) = DerivativePolyCoeffient_helper(n_order, 0, 0);
        beq_wp(mid_idx-1) = waypoints(mid_idx);
    end

    % carry_idx = carry_idx+1;
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    %
    %
    %
    %
    for mid_idx = 1 : n_seg-1
        Aeq_con_p(mid_idx, (mid_idx-1)*8+1:(mid_idx)*8) = DerivativePolyCoeffient_helper(n_order, 0, ts(mid_idx));
        Aeq_con_p(mid_idx, (mid_idx)*8+1:(mid_idx+1)*8) = -DerivativePolyCoeffient_helper(n_order, 0, 0);
    end
    
    % carry_idx = carry_idx +2;

    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    %
    %
    %
    %
    for mid_idx = 1 : n_seg-1
        Aeq_con_v(mid_idx, (mid_idx-1)*8+1:(mid_idx)*8) = DerivativePolyCoeffient_helper(n_order, 1, ts(mid_idx));
        Aeq_con_v(mid_idx, (mid_idx)*8+1:(mid_idx+1)*8) = -DerivativePolyCoeffient_helper(n_order, 1, 0);
    end

    % carry_idx = carry_idx +2;

    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    %
    %
    %
    %
    for mid_idx = 1 : n_seg-1
        Aeq_con_a(mid_idx, (mid_idx-1)*8+1:(mid_idx)*8) = DerivativePolyCoeffient_helper(n_order, 2, ts(mid_idx));
        Aeq_con_a(mid_idx, (mid_idx)*8+1:(mid_idx+1)*8) = -DerivativePolyCoeffient_helper(n_order, 2, 0);
    end

    % carry_idx = carry_idx +2;

    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    %
    %
    %
    %
    for mid_idx = 1 : n_seg-1
        Aeq_con_j(mid_idx, (mid_idx-1)*8+1:(mid_idx)*8) = DerivativePolyCoeffient_helper(n_order, 3, ts(mid_idx));
        Aeq_con_j(mid_idx, (mid_idx)*8+1:(mid_idx+1)*8) = -DerivativePolyCoeffient_helper(n_order, 3, 0);
    end


    size(Aeq_con_p)
    size(Aeq_con_v)
    size(Aeq_con_a)
    size(Aeq_con_j)
    size(Aeq_start)
    
    size(Aeq_end)
    size(Aeq_wp)
    
    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    size(Aeq_con)
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    a = 0.8 
    size(beq_start)
    size(beq_end)
    size(beq_wp)
    size(beq_con)
    beq = [beq_start; beq_end; beq_wp; beq_con];
end