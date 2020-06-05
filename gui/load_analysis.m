function handles = load_analysis(handles, varargin)

% Loop through input params
nVarIn = length(varargin);

for i = 1:nVarIn

    if strcmpi(varargin{i}, 'Subjects')
        
        handles.subjects = varargin{i+1};
    
    elseif strcmpi(varargin{i}, 'Params')
        
        handles.params = varargin{i+1};
    
    elseif strcmpi(varargin{i}, 'Analysis')
        
        handles.GDdims = varargin{i+3};
        handles = load_GDdims(handles, varargin{i+1}, handles.NM.label, handles.GDdims);
        
    elseif strcmpi(varargin{i}, 'Visdata')
    
        vis = varargin{i+1};
        if ~isempty( vis  ) 
            if numel( vis )>1
                for n=1:numel( vis )
                    try
                        handles.visdata{n} = vis{n}{1};
                    catch
                        handles.visdata{n} = vis{n};
                    end
                    handles.visdata_table(n) = create_visdata_tables(handles.visdata{n}, [], [], 'create');
                end
            else
                for n=1:numel( vis{1} )
                    try
                        handles.visdata{n} = vis{1}{n};
                    catch
                        handles.visdata{n} = vis{n};
                    end
                    handles.visdata_table(n) = create_visdata_tables(handles.visdata{n}, [], [], 'create');
                end
            end
        end
        
    elseif strcmpi(varargin{i}, 'OOCVdata')
        
        handles = load_OOCV(handles, varargin{i+1});
    
    end
    
end
