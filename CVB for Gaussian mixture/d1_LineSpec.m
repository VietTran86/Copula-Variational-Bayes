%////////////////////////////////////////////////////////////////////////////////////////////////
st.Line  = 'LineWidth';
st.Edge  = 'MarkerEdgeColor';
st.Face  = 'MarkerFaceColor';
st.Msize = 'MarkerSize';
st.Disp  = 'DisplayName'; 
%--------------------------------------------
spec.kmean = {'k',   st.Disp,'k-means'          st.Line,1,st.Msize,10};
spec.emL   = {'xb:',   st.Disp,'EM_1'    st.Line,1,st.Msize,10};
spec.emMU  = {'ob:',   st.Disp,'EM_2'     st.Line,1,st.Msize,10};
spec.VB    = {'db--',   st.Disp,'VB'       st.Line,1,st.Msize,10};
spec.CVB1  = {'hc-.',    st.Disp,'CVB_1'    st.Line,1,st.Msize,10};
spec.CVB2  = {'sm-.',    st.Disp,'CVB_2'      st.Line,1,st.Msize,10};
spec.CVB3  = {'*r-.',    st.Disp,'CVB_3'       st.Line,1,st.Msize,10};
%////////////////////////////////////////////////////////////////////////////////////////////////