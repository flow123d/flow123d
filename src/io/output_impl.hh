#ifndef OUTPUT_IMPL_HH
#define OUTPUT_IMPL_HH




template<int spacedim, class Value>
void OutputTime::register_data(const Input::Record &in_rec,
        const RefType type,
        MultiField<spacedim, Value> &multi_field)
{
    for (unsigned long index=0; index < multi_field.size(); index++) {
        OutputTime::register_data(in_rec, type, multi_field[index] );
    }
}


template<int spacedim, class Value>
void OutputTime::register_data(const Input::Record &in_rec,
        const RefType ref_type,
        Field<spacedim, Value> &field_ref)
{
	Field<spacedim,Value> *field = &field_ref;
    string name_ = field->name();
    OutputDataBase *output_data;
    unsigned int item_count = 0, comp_count = 0, node_id;

    // TODO: do not try to find empty string and raise exception

    // Try to find record with output stream (the key is name of data)
    Input::Iterator<string> stream_name_iter = in_rec.find<string>(name_);

    // If record was not found, then throw exception
    if(!stream_name_iter) {
        return;
    }

    // Try to find existing output stream
    OutputTime *output_time = OutputTime::output_stream_by_name(*stream_name_iter);

    /* It's possible now to do output to the file only in the first process */
    if(output_time == NULL || output_time->rank != 0) {
        /* TODO: do something, when support for Parallel VTK is added */
        return;
    }

    // TODO: remove const_cast after resolving problems with const Mesh.
    Mesh *mesh = const_cast<Mesh *>(field->mesh());

    if(output_time->get_mesh() == NULL) {
        output_time->set_mesh(mesh);
    }

    ElementFullIter ele = ELEMENT_FULL_ITER(mesh, NULL);
    Node *node;
    int corner_index = 0;
    int node_index = 0;
    int ele_index = 0;

    output_data = output_time->output_data_by_field((FieldCommon*)field,
            ref_type);

    switch(ref_type) {
    case NODE_DATA:
        item_count = mesh->n_nodes();
        break;
    case CORNER_DATA:
        // Compute number of all corners
        item_count = 0;
        FOR_ELEMENTS(mesh, ele) {
            item_count += ele->n_nodes();
        }
        break;
    case ELEM_DATA:
        item_count = mesh->n_elements();
        break;
    }

    if(output_data == NULL) {
        /* This is problematic part, because of templates :-( */
        if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
            output_data = (OutputDataBase*)new OutputData<int>(field, item_count, 1);
        } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
            output_data = (OutputDataBase*)new OutputData<int>(field, item_count, 3);
        } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
            output_data = (OutputDataBase*)new OutputData<unsigned int>(field, item_count, 1);
        } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
            output_data = (OutputDataBase*)new OutputData<unsigned int>(field, item_count, 3);
        } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
            output_data = (OutputDataBase*)new OutputData<double>(field, item_count, 1);
        } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
            output_data = (OutputDataBase*)new OutputData<double>(field, item_count, 3);
        } else {
            /* TODO: this will not be necessary */
            throw "Try to register unsupported data type.";
        }

        switch(ref_type) {
        case NODE_DATA:
            output_time->node_data.push_back(output_data);
            break;
        case CORNER_DATA:
            output_time->corner_data.push_back(output_data);
            break;
        case ELEM_DATA:
            output_time->elem_data.push_back(output_data);
            break;
        }
    }

    unsigned int *count = new unsigned int[item_count];

    /* Copy data to array */
    switch(ref_type) {
    case NODE_DATA:

        // Initialize arrays
        for(node_id=0; node_id<item_count; node_id++) {
            if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                (*(OutputData<int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                (*(OutputData<int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                (*(OutputData<double>*)output_data)[node_id] = 0;
            } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                (*(OutputData<double>*)output_data)[node_id] = 0;
            }
            count[node_id] = 0;
        }

        /* Copy data to temporary array */
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, node_id) {
                node = ele->node[node_id];
                node_index = mesh->node_vector.index(ele->node[node_id]);
                if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                    (*(OutputData<int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                    (*(OutputData<int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                    (*(OutputData<unsigned int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                    (*(OutputData<unsigned int>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                    (*(OutputData<double>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                    (*(OutputData<double>*)output_data)[node_index] += field->value(node->point(), mesh->element_accessor(ele_index));
                }
                count[mesh->node_vector.index(ele->node[node_id])]++;
            }
        }

        // Compute mean values at nodes
        for(node_id=0; node_id<item_count; node_id++) {
            if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                (*(OutputData<int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                (*(OutputData<int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                (*(OutputData<unsigned int>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                (*(OutputData<double>*)output_data)[node_id] /= count[node_id];
            } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                (*(OutputData<double>*)output_data)[node_id] /= count[node_id];
            }
        }

        break;
    case CORNER_DATA:
        FOR_ELEMENTS(mesh, ele) {
            FOR_ELEMENT_NODES(ele, node_id) {
                node = ele->node[node_id];
                if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                    (*(OutputData<int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                    (*(OutputData<int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                    (*(OutputData<unsigned int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                    (*(OutputData<unsigned int>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                    (*(OutputData<double>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                    (*(OutputData<double>*)output_data)[corner_index] = field->value(node->point(), mesh->element_accessor(ele_index));
                }
                corner_index++;
            }
            ele_index++;
        }
        break;
    case ELEM_DATA:
        FOR_ELEMENTS(mesh, ele) {
            if(typeid(Value) == typeid(FieldValue<1>::Integer)) {
                (*(OutputData<int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::IntVector)) {
                (*(OutputData<int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::Enum)) {
                (*(OutputData<unsigned int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::EnumVector)) {
                (*(OutputData<unsigned int>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::Scalar)) {
                (*(OutputData<double>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            } else if(typeid(Value) == typeid(FieldValue<1>::Vector)) {
                (*(OutputData<double>*)output_data)[ele_index] = field->value(ele->centre(), mesh->element_accessor(ele_index));
            }
            ele_index++;
        }
        break;
    }

    /* Set the last time */
    if(output_time->time < field->time()) {
        output_time->time = field->time();
    }
}


#endif
