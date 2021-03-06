//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Alexandru Ichim
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#include "viewer.h"

void Viewer::select_point(const Eigen::Vector2i & pixel) {
	Eigen::Matrix4f model, view, projection;
	computeCameraMatrices(model, view, projection);
	Matrix4f MVP = projection * view * model;
	Matrix4f invMVP = MVP.inverse();

	Eigen::Vector3f origin = nanogui::unproject(Eigen::Vector3f(pixel[0], pixel[1], 0), view * model, projection, mSize);
	Eigen::Vector3f endpoint = nanogui::unproject(Eigen::Vector3f(pixel[0], pixel[1], 1), view * model, projection, mSize);
	Eigen::Vector3f direction = (endpoint - origin) / (endpoint - origin).norm();

	Eigen::Vector3f closest_vertex = mesh_->get_closest_vertex(origin, direction, contraint_indices_[edit_constraint_index_]);
	mesh_->set_selection(closest_vertex, edit_constraint_index_);
	cout << contraint_indices_[0] << " " << contraint_indices_[1] << " " << contraint_indices_[2] << " " << contraint_indices_[3] << endl;
}

bool Viewer::keyboardEvent(int key, int scancode, int action, int modifiers) {
	if (Screen::keyboardEvent(key, scancode, action, modifiers)) {
		return true;
	}
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		setVisible(false);
		return true;
	}
	return false;
}

void Viewer::draw(NVGcontext *ctx) {
	/* Draw the user interface */
	Screen::draw(ctx);
}

Vector2f Viewer::getScreenCoord() {
	Vector2i pos = mousePos();
	return Vector2f(2.0f * (float)pos.x() / width() - 1.0f,
		1.0f - 2.0f * (float)pos.y() / height());
}

void Viewer::drawContents() {
	using namespace nanogui;

	/* Draw the window contents using OpenGL */
	shader_.bind();

	Eigen::Matrix4f model, view, proj;
	computeCameraMatrices(model, view, proj);

	Matrix4f mv = view*model;
	Matrix4f p = proj;

	/* MVP uniforms */
	shader_.setUniform("MV", mv);
	shader_.setUniform("P", p);

	// Setup OpenGL (making sure the GUI doesn't disable these
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	// Render everything
	if (wireframe_) {
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	Vector3f colors(1.0, 1.0, 1.0);
	shader_.setUniform("intensity", colors);
	if (color_mode == CURVATURE) {
		shader_.setUniform("color_mode", int(curvature_type));
	}
	else {
		shader_.setUniform("color_mode", int(color_mode));
	}
	shader_.drawIndexed(GL_TRIANGLES, 0, mesh_->get_number_of_face());

	if (wireframe_) {
		glDisable(GL_POLYGON_OFFSET_FILL);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		colors << 0.6, 0.6, 0.6;
		shader_.setUniform("intensity", colors);
		shader_.drawIndexed(GL_TRIANGLES, 0, mesh_->get_number_of_face());
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	if (normals_) {
		shaderNormals_.bind();
		shaderNormals_.setUniform("MV", mv);
		shaderNormals_.setUniform("P", p);
		shaderNormals_.drawIndexed(GL_TRIANGLES, 0, mesh_->get_number_of_face());
	}

	if (isolines_) {
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(2.0);
		shaderIsolines_.bind();
		shaderIsolines_.setUniform("MV", mv);
		shaderIsolines_.setUniform("P", p);
		shaderIsolines_.drawArray(GL_LINES, 0, mesh_->isolines_points_.size());
		glDisable(GL_LINE_SMOOTH);
	}

	shaderSelection_.bind();
	shaderSelection_.setUniform("MV", mv);
	shaderSelection_.setUniform("P", p);
	glEnable(GL_PROGRAM_POINT_SIZE);
	shaderSelection_.drawIndexed(GL_POINTS, 0, mesh_->get_number_of_face());
	glDisable(GL_PROGRAM_POINT_SIZE);		
	
}

bool Viewer::scrollEvent(const Vector2i &p, const Vector2f &rel) {
	if (!Screen::scrollEvent(p, rel)) {
		camera_.zoom = max(0.1, camera_.zoom * (rel.y() > 0 ? 1.1 : 0.9));
	}
	return true;
}

bool Viewer::mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
	int button, int modifiers) {

	if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
		if (camera_.arcball.motion(p)) {
			//
		}
		else if (translate_) {
			Eigen::Matrix4f model, view, proj;
			computeCameraMatrices(model, view, proj);
			Point mesh_center = mesh_->get_mesh_center();
			float zval = nanogui::project(Vector3f(mesh_center.x,
				mesh_center.y,
				mesh_center.z),
				view * model, proj, mSize).z();
			Eigen::Vector3f pos1 = nanogui::unproject(
				Eigen::Vector3f(p.x(), mSize.y() - p.y(), zval),
				view * model, proj, mSize);
			Eigen::Vector3f pos0 = nanogui::unproject(
				Eigen::Vector3f(translateStart_.x(), mSize.y() -
					translateStart_.y(), zval), view * model, proj, mSize);
			camera_.modelTranslation = camera_.modelTranslation_start + (pos1 - pos0);
		}
	}
	return true;
}

bool Viewer::mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) {
	if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
		if (button == GLFW_MOUSE_BUTTON_1 && modifiers == 0) {
			camera_.arcball.button(p, down);
		}
		else if (button == GLFW_MOUSE_BUTTON_2 ||
			(button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_SHIFT)) {
			camera_.modelTranslation_start = camera_.modelTranslation;
			translate_ = true;
			translateStart_ = p;
		}
		else if (button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_CONTROL) {
			select_point(Eigen::Vector2i(p.x(), mSize.y() - p.y()));
			refresh_selection();
		}
	}
	if (button == GLFW_MOUSE_BUTTON_1 && !down) {
		camera_.arcball.button(p, false);
	}
	if (!down) {
		translate_ = false;
	}
	return true;
}

void Viewer::initShaders() {
	// Shaders
	shader_.init(
		"a_simple_shader",

		/* Vertex shader */
		"#version 330\n"
		"uniform mat4 MV;\n"
		"uniform mat4 P;\n"
		"uniform int color_mode;\n"
		"uniform vec3 intensity;\n"

		"in vec3 position;\n"
		"in vec3 valence_color;\n"
		"in vec3 unicruvature_color;\n"
		"in vec3 laplacian_color;\n"
		"in vec3 curvature_color;\n"
		"in vec3 gaussian_curv_color;\n"
		"in vec3 normal;\n"

		"out vec3 fcolor;\n"
		"out vec3 fnormal;\n"
		"out vec3 view_dir;\n"
		"out vec3 light_dir;\n"

		"void main() {\n"
		"    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
		"    gl_Position = P * vpoint_mv;\n"
		"    if (color_mode == 1) {\n"
		"        fcolor = valence_color;\n"
		"    } else if (color_mode == 2) {\n"
		"        fcolor = unicruvature_color;\n"
		"    } else if (color_mode == 3) {\n"
		"        fcolor = curvature_color;\n"
		"    } else if (color_mode == 4) {\n"
		"        fcolor = gaussian_curv_color;\n"
		"    } else if (color_mode == 5) {\n"
		"        fcolor = laplacian_color;\n"
		"    } else {\n"
		"        fcolor = intensity;\n"
		"    }\n"
		"    fnormal = mat3(transpose(inverse(MV))) * normal;\n"
		"    light_dir = vec3(0.0, 3.0, 3.0) - vpoint_mv.xyz;\n"
		"    view_dir = -vpoint_mv.xyz;\n"
		"}",

		/* Fragment shader */
		"#version 330\n"
		"uniform int color_mode;\n"
		"uniform vec3 intensity;\n"

		"in vec3 fcolor;\n"
		"in vec3 fnormal;\n"
		"in vec3 view_dir;\n"
		"in vec3 light_dir;\n"

		"out vec4 color;\n"

		"void main() {\n"
		"    vec3 c = vec3(0.0);\n"
		"    if (color_mode >= 0) {\n"
		"        c += vec3(1.0)*vec3(0.15, 0.15, 0.15);\n"
		"        vec3 n = normalize(fnormal);\n"
		"        vec3 v = normalize(view_dir);\n"
		"        vec3 l = normalize(light_dir);\n"
		"        float lambert = dot(n,l);\n"
		"        if(lambert > 0.0) {\n"
		"            c += vec3(1.0)*vec3(0.85, 0.85, 0.85)*lambert;\n"
		"            vec3 v = normalize(view_dir);\n"
		"            vec3 r = reflect(-l,n);\n"
		"            c += vec3(1.0)*vec3(0.7, 0.7, 0.7)*pow(max(dot(r,v), 0.0), 90.0);\n"
		"        }\n"
		"        c *= fcolor;\n"
		"    } else {\n"
		"       c = fcolor;\n"
		"    }\n"
		"    if (intensity == vec3(0.0)) {\n"
		"        c = intensity;\n"
		"    }\n"
		"    color = vec4(c, 1.0);\n"
		"}"
	);

	shaderNormals_.init(
		"normal_shader",
		/* Vertex shader */
		"#version 330\n\n"
		"in vec3 position;\n"
		"in vec3 normal;\n"
		"uniform mat4 MV;\n"
		"uniform mat4 P;\n"
		"uniform int normal_selector;\n"
		"out VS_OUT {\n"
		"    mat3 normal_mat;\n"
		"    vec3 normal;\n"
		"} vs_out;\n"
		"void main() {\n"
		"  gl_Position = vec4(position, 1.0);\n"
		"    vs_out.normal = normal;\n"
		"    vs_out.normal_mat = mat3(transpose(inverse(MV)));\n"
		"}",
		/* Fragment shader */
		"#version 330\n\n"
		"out vec4 frag_color;\n"
		"void main() {\n"
		"   frag_color = vec4(0.0, 1.0, 0.0, 1.0);\n"
		"}",
		/* Geometry shader */
		"#version 330\n\n"
		"layout (triangles) in;\n"
		"layout (line_strip, max_vertices = 6) out;\n"
		"uniform mat4 MV;\n"
		"uniform mat4 P;\n"
		"in VS_OUT {\n"
		"    mat3 normal_mat;\n"
		"    vec3 normal;\n"
		"} gs_in[];\n"
		"void createline(int index) {\n"
		"   gl_Position = P * MV * gl_in[index].gl_Position;\n"
		"   EmitVertex();\n"
		"   vec4 normal_mv = vec4(normalize(gs_in[index].normal_mat *\n"
		"                                   gs_in[index].normal), 1.0f);\n"
		"   gl_Position = P * (MV * gl_in[index].gl_Position\n"
		"                      + normal_mv * 0.035f);\n"
		"   EmitVertex();\n"
		"   EndPrimitive();\n"
		"}\n"
		"void main() {\n"
		"   createline(0);\n"
		"   createline(1);\n"
		"   createline(2);\n"
		"}"
	);

	shaderSelection_.init(
	"selection_shader",

	"#version 330\n"
	"in vec3 position;\n"
	"uniform mat4 MV;\n"
	"uniform mat4 P;\n"
	"void main() {\n"
	"    vec4 vpoint_mv = MV * vec4(position, 1.0);\n"
	"    gl_Position = P * vpoint_mv;\n"
	"    gl_PointSize = 10.0;\n"
	"}",

	"#version 330\n"
	"out vec4 color;\n"
	"void main() {\n"
	"	 if (gl_PrimitiveID == 0) color = vec4(0.7, 0.0, 0.2, 1.0);\n"
	"	 if (gl_PrimitiveID == 1) color = vec4(0.3, 0.2, 0.7, 1.0);\n"
	"	 if (gl_PrimitiveID == 2) color = vec4(0.13, 0.7, 0.3, 1.0);\n"
	"	 if (gl_PrimitiveID == 3) color = vec4(1, 0.5, 0.15, 1.0);\n"
	"}"
	);

	shaderIsolines_.init(
		"isolines_shader",
		/* Vertex shader */

		"#version 330\n\n"
		"in vec3 position;\n"
		"uniform mat4 MV;\n"
		"uniform mat4 P;\n"
		"void main() {\n"
		"  gl_Position = P*MV*vec4(position, 1.0);\n"
		"}",
		/* Fragment shader */
		"#version 330\n\n"
		"out vec4 color;\n"
		"void main() {\n"
		"   color = vec4(1.0, 1.0, 1.0, 0.0);\n"
		"}"
	);
}

Viewer::Viewer() : nanogui::Screen(Eigen::Vector2i(1024, 768), "DGP Viewer") {

	window_ = new Window(this, "Controls");
	window_->setPosition(Vector2i(15, 15));
	window_->setLayout(new GroupLayout());

	PopupButton *popupBtn = new PopupButton(window_, "Open a mesh", ENTYPO_ICON_EXPORT);
	Popup *popup = popupBtn->popup();
	popup->setLayout(new GroupLayout());

	Button* b = new Button(popup, "Bunny");
	b->setCallback([this]() {
		mesh_->load_mesh("../data/bunny.off");
		this->refresh_mesh();
		this->refresh_trackball_center();
	});
	b = new Button(popup, "Max-Planck");
	b->setCallback([this]() {
		mesh_->load_mesh("../data/max.off");
		this->refresh_mesh();
		this->refresh_trackball_center();
	});

	b = new Button(popup, "Open mesh ...");
	b->setCallback([this]() {
		string filename = nanogui::file_dialog({ { "obj", "Wavefront OBJ" },
		{ "ply", "Stanford PLY" },
		{ "aln", "Aligned point cloud" },
		{ "off", "Object File Format" }
		}, false);
		if (filename != "") {
			mesh_->load_mesh(filename);
			this->refresh_mesh();
			this->refresh_trackball_center();
		}
	});

	new Label(window_, "Display Control", "sans-bold");

	b = new Button(window_, "Wireframe");
	b->setFlags(Button::ToggleButton);
	b->setChangeCallback([this](bool wireframe) {
		this->wireframe_ = !this->wireframe_;
	});
	b = new Button(window_, "Normals");
	b->setFlags(Button::ToggleButton);
	b->setChangeCallback([this](bool normals) {
		this->normals_ = !this->normals_;
	});

	b = new Button(window_, "Laplacian");
	b->setFlags(Button::ToggleButton);
	b->setChangeCallback([this](bool laplacian) {
		if (laplacian) {
			this->color_mode = LAPLACIAN;
		}
		else {
			this->color_mode = NORMAL;
		}
	});


	b = new Button(window_, "Isolines");
	b->setFlags(Button::ToggleButton);
	b->setChangeCallback([this](bool iso_curve) {
		this->isolines_ = !this->isolines_;
	});

	new Label(window_, "Constraints", "sans-bold");
	popupBtn = new PopupButton(window_, "Edit Constraint");
	popup = popupBtn->popup();
	popup->setLayout(new GroupLayout());
	b = new Button(popup, "Number 0");
	b->setCallback([this]() {
		edit_constraint_index_ = 0;
	});
	b = new Button(popup, "Number 1");
	b->setCallback([this]() {
		edit_constraint_index_ = 1;
	});
	b = new Button(popup, "Number 2");
	b->setCallback([this]() {
		edit_constraint_index_ = 2;
	});
	b = new Button(popup, "Number 3");
	b->setCallback([this]() {
		edit_constraint_index_ = 3;
	});

	b = new Button(window_, "Harmonic Functions");
	b->setChangeCallback([this](bool normals) {
		std::vector<size_t> constrain_indices_first = { contraint_indices_[0], contraint_indices_[1] };
	    std::vector<size_t> constrain_indices_second = { contraint_indices_[2], contraint_indices_[3] };
		mesh_->harmonic_function(constrain_indices_first, "v:harmonic_function_0");
		mesh_->harmonic_function(constrain_indices_second, "v:harmonic_function_1");
		cout << "DONE" << endl;
		mesh_->compute_mesh_properties();
		this->refresh_mesh();
	});

	new Label(window_, "Number of intervals:", "sans-bold");
	intervalsTextBox = new IntBox<int>(window_, 10);
	intervalsTextBox->setEditable(true);
	intervalsTextBox->setFixedSize(Vector2i(50, 20));
	intervalsTextBox->setDefaultValue("10");
	intervalsTextBox->setFontSize(16);

	b = new Button(window_, "Compute Isolines");
	b->setChangeCallback([this](bool normals) {
		std::vector<size_t> constrain_indices_first = { contraint_indices_[0], contraint_indices_[1] };
		std::vector<size_t> constrain_indices_second = { contraint_indices_[2], contraint_indices_[3] };
		mesh_->isolines_points_.clear();
		mesh_->compute_isolines(constrain_indices_first, "v:harmonic_function_0", intervalsTextBox->value());
		mesh_->compute_isolines(constrain_indices_second, "v:harmonic_function_1", intervalsTextBox->value());
		cout << "DONE" << endl;
		mesh_->compute_mesh_properties();
		this->refresh_mesh();
	});

	performLayout();

	initShaders();
	mesh_ = new mesh_processing::MeshProcessing("../data/bunny_1000.obj");
	this->refresh_mesh();
	this->refresh_trackball_center();
}

void Viewer::refresh_trackball_center() {
	// Re-center the mesh
	Point mesh_center = mesh_->get_mesh_center();
	camera_.arcball = Arcball();
	camera_.arcball.setSize(mSize);
	camera_.modelZoom = 2 / mesh_->get_dist_max();
	camera_.modelTranslation = -Vector3f(mesh_center.x, mesh_center.y, mesh_center.z);
}

void Viewer::refresh_mesh() {
	shader_.bind();
	shader_.uploadIndices(*(mesh_->get_indices()));
	shader_.uploadAttrib("position", *(mesh_->get_points()));
	shader_.uploadAttrib("valence_color", *(mesh_->get_colors_valence()));
	shader_.uploadAttrib("unicruvature_color", *(mesh_->get_colors_unicurvature()));
	shader_.uploadAttrib("curvature_color", *(mesh_->get_color_curvature()));
	shader_.uploadAttrib("gaussian_curv_color", *(mesh_->get_colors_gaussian_curv()));
	shader_.uploadAttrib("laplacian_color", *(mesh_->get_colors_laplacian()));
	shader_.uploadAttrib("normal", *(mesh_->get_normals()));
	shader_.setUniform("color_mode", int(color_mode));
	shader_.setUniform("intensity", Vector3f(1.0, 1.0, 1.0));

	shaderNormals_.bind();
	shaderNormals_.shareAttrib(shader_, "indices");
	shaderNormals_.shareAttrib(shader_, "position");
	shaderNormals_.shareAttrib(shader_, "normal");

	refresh_selection();

	refresh_isolines();
}

void Viewer::refresh_selection() {
	shaderSelection_.bind();
	Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> indices = Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic>(4, 1);
	indices << 0, 1, 2, 3;
	Eigen::MatrixXf selection = (*(mesh_->get_selection()));
	shaderSelection_.uploadIndices(indices);
	shaderSelection_.uploadAttrib("position", selection);
}

void Viewer::refresh_isolines() {
	shaderIsolines_.bind();
	MatrixXf isolines_matrix(3, mesh_->isolines_points_.size());
	for (size_t j = 0; j < mesh_->isolines_points_.size(); j++)
		isolines_matrix.col(j) <<
		mesh_->isolines_points_[j].x,
		mesh_->isolines_points_[j].y,
		mesh_->isolines_points_[j].z;

	shaderIsolines_.uploadAttrib("position", isolines_matrix);
}

void Viewer::computeCameraMatrices(Eigen::Matrix4f &model,
	Eigen::Matrix4f &view,
	Eigen::Matrix4f &proj) {

	view = nanogui::lookAt(camera_.eye, camera_.center, camera_.up);

	float fH = std::tan(camera_.viewAngle / 360.0f * M_PI) * camera_.dnear;
	float fW = fH * (float)mSize.x() / (float)mSize.y();

	proj = nanogui::frustum(-fW, fW, -fH, fH, camera_.dnear, camera_.dfar);
	model = camera_.arcball.matrix();

	model = nanogui::scale(model, Eigen::Vector3f::Constant(camera_.zoom * camera_.modelZoom));
	model = nanogui::translate(model, camera_.modelTranslation);
}

Viewer::~Viewer() {
	shader_.free();
	shaderNormals_.free();
}
